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
#ifndef TESTPOTTSUPDATERULES_HPP_
#define TESTPOTTSUPDATERULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractPottsUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CellwiseData.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "PottsMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"

class TestPottsUpdateRules : public AbstractCellBasedTestSuite
{
public:

    void TestVolumeConstraintPottsUpdateRuleMethods() throw (Exception)
    {
		// Create a simple 2D PottsMesh with 2 elements
		PottsMeshGenerator<2> generator(4, 1, 2, 4, 2, 2, 1, 1, 1, true); // last bool makes elements start in bottom left
		PottsMesh<2>* p_mesh = generator.GetMesh();

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

		// Create cell population
		PottsBasedCellPopulation cell_population(*p_mesh, cells);

		// Create an update law system
		VolumeConstraintPottsUpdateRule<2> volume_constraint;

		// Test get/set methods
		TS_ASSERT_DELTA(volume_constraint.GetDeformationEnergyParameter(), 0.5, 1e-12);
		TS_ASSERT_DELTA(volume_constraint.GetMatureCellTargetVolume(), 16.0, 1e-12);

		volume_constraint.SetDeformationEnergyParameter(1.0);
		volume_constraint.SetMatureCellTargetVolume(0.6);

		TS_ASSERT_DELTA(volume_constraint.GetDeformationEnergyParameter(), 1.0, 1e-12);
		TS_ASSERT_DELTA(volume_constraint.GetMatureCellTargetVolume(), 0.6, 1e-12);

		volume_constraint.SetDeformationEnergyParameter(0.5);
		volume_constraint.SetMatureCellTargetVolume(4.0);

		// Test EvaluateHamiltonianContribution()

		double alpha = volume_constraint.GetDeformationEnergyParameter();

		// Both points lie within cell 0
		TS_ASSERT_THROWS_THIS(volume_constraint.EvaluateHamiltonianContribution(0, 1, cell_population),
		                      "The current node and target node must not be in the same element.");

		// Both points lie within cell 1
		TS_ASSERT_THROWS_THIS(volume_constraint.EvaluateHamiltonianContribution(8, 9, cell_population),
							  "The current node and target node must not be in the same element.");

		// Both points lie within cell medium
		TS_ASSERT_THROWS_THIS(volume_constraint.EvaluateHamiltonianContribution(2, 3, cell_population),
		                      "At least one of the current node or target node must be in an element.");

		// Current site in cell 0; target site in cell medium
		double contribution = volume_constraint.EvaluateHamiltonianContribution(5, 6, cell_population);
		TS_ASSERT_DELTA(contribution, alpha, 1e-6);

		// Current site in cell medium; target site in cell 0
		contribution = volume_constraint.EvaluateHamiltonianContribution(6, 5, cell_population);
		TS_ASSERT_DELTA(contribution, alpha, 1e-6);

		// Current site in cell 0; target site in cell 1
		contribution = volume_constraint.EvaluateHamiltonianContribution(5, 9, cell_population);
		TS_ASSERT_DELTA(contribution, 2.0*alpha, 1e-6);

	}

    void TestArchiveVolumeConstraintPottsUpdateRule() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "VolumeConstraintPottsUpdateRule.arch";

        {
            VolumeConstraintPottsUpdateRule<2> update_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            update_rule.SetDeformationEnergyParameter(0.5);
            update_rule.SetMatureCellTargetVolume(0.6);

            // Serialize via pointer to most abstract class possible
            AbstractPottsUpdateRule<2>* const p_update_rule = &update_rule;
            output_arch << p_update_rule;
        }

        {
            AbstractPottsUpdateRule<2>* p_update_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_update_rule;

            // Test the member data
            TS_ASSERT_DELTA((static_cast<VolumeConstraintPottsUpdateRule<2>*>(p_update_rule))->GetDeformationEnergyParameter(), 0.5, 1e-6);
            TS_ASSERT_DELTA((static_cast<VolumeConstraintPottsUpdateRule<2>*>(p_update_rule))->GetMatureCellTargetVolume(), 0.6, 1e-6);

            // Tidy up
            delete p_update_rule;
        }
    }
    
    void TestSurfaceAreaConstraintPottsUpdateRuleMethods() throw (Exception)
    {
        // Create a simple 2D PottsMesh with 2 elements
        PottsMeshGenerator<2> generator(4, 1, 2, 4, 2, 2, 1, 1, 1, true); // last bool makes elements start in bottom left
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation cell_population(*p_mesh, cells);

        // Create an update law system
        SurfaceAreaConstraintPottsUpdateRule<2> surface_area_constraint;

        // Test get/set methods
        TS_ASSERT_DELTA(surface_area_constraint.GetDeformationEnergyParameter(), 0.5, 1e-12);
        TS_ASSERT_DELTA(surface_area_constraint.GetMatureCellTargetSurfaceArea(), 16.0, 1e-12);

        surface_area_constraint.SetDeformationEnergyParameter(1.0);
        surface_area_constraint.SetMatureCellTargetSurfaceArea(0.6);

        TS_ASSERT_DELTA(surface_area_constraint.GetDeformationEnergyParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(surface_area_constraint.GetMatureCellTargetSurfaceArea(), 0.6, 1e-12);

        surface_area_constraint.SetDeformationEnergyParameter(0.5);
        surface_area_constraint.SetMatureCellTargetSurfaceArea(8.0);

		// Test EvaluateHamiltonianContribution()
        double gamma = surface_area_constraint.GetDeformationEnergyParameter();

		// Both points lie within cell 0
		TS_ASSERT_THROWS_THIS(surface_area_constraint.EvaluateHamiltonianContribution(0, 1, cell_population),
		                      "The current node and target node must not be in the same element.");

		// Both points lie within cell 1
		TS_ASSERT_THROWS_THIS(surface_area_constraint.EvaluateHamiltonianContribution(8, 9, cell_population),
							  "The current node and target node must not be in the same element.");

		// Both points lie within cell medium
		TS_ASSERT_THROWS_THIS(surface_area_constraint.EvaluateHamiltonianContribution(2, 3, cell_population),
		                       "At least one of the current node or target node must be in an element.");

		// Current site in cell 0; target site in cell medium Cell on edge of domain
		double contribution = surface_area_constraint.EvaluateHamiltonianContribution(1, 2, cell_population);
		TS_ASSERT_DELTA(contribution, 4.0*gamma, 1e-6);

		// Current site in cell 0; target site in cell medium
		contribution = surface_area_constraint.EvaluateHamiltonianContribution(5, 6, cell_population);
		TS_ASSERT_DELTA(contribution, 4.0*gamma, 1e-6);

		// Current site in cell medium; target site in cell 0
		contribution = surface_area_constraint.EvaluateHamiltonianContribution(6, 5, cell_population);
		TS_ASSERT_DELTA(contribution, 0.0, 1e-6);

		// Current site in cell 0; target site in cell 1
		contribution = surface_area_constraint.EvaluateHamiltonianContribution(5, 9, cell_population);
		TS_ASSERT_DELTA(contribution, 4.0*gamma, 1e-6);
    }

    void TestArchiveSurfaceAreaConstraintPottsUpdateRule() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SurfaceAreaConstraintPottsUpdateRule.arch";

        {
            SurfaceAreaConstraintPottsUpdateRule<2> update_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            update_rule.SetDeformationEnergyParameter(0.5);
            update_rule.SetMatureCellTargetSurfaceArea(0.6);

            // Serialize via pointer to most abstract class possible
            AbstractPottsUpdateRule<2>* const p_update_rule = &update_rule;
            output_arch << p_update_rule;
        }

        {
            AbstractPottsUpdateRule<2>* p_update_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_update_rule;

            // Test the member data
            TS_ASSERT_DELTA((static_cast<SurfaceAreaConstraintPottsUpdateRule<2>*>(p_update_rule))->GetDeformationEnergyParameter(), 0.5, 1e-6);
            TS_ASSERT_DELTA((static_cast<SurfaceAreaConstraintPottsUpdateRule<2>*>(p_update_rule))->GetMatureCellTargetSurfaceArea(), 0.6, 1e-6);

            // Tidy up
            delete p_update_rule;
        }
    }

    void TestAdhesionPottsUpdateRuleMethods() throw (Exception)
    {
    	// Create a simple 2D PottsMesh with 2 elements
    	PottsMeshGenerator<2> generator(4, 1, 2, 4, 2, 2, 1, 1, 1, true); // last bool makes elements start in bottom left
		PottsMesh<2>* p_mesh = generator.GetMesh();

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

		// Create cell population
		PottsBasedCellPopulation cell_population(*p_mesh, cells);

		// Create an update law system
		AdhesionPottsUpdateRule<2> adhesion_update;

		// Test get/set methods
	 	TS_ASSERT_DELTA(adhesion_update.GetCellCellAdhesionEnergyParameter(), 0.1, 1e-12);
		TS_ASSERT_DELTA(adhesion_update.GetCellBoundaryAdhesionEnergyParameter(), 0.2, 1e-12);

		adhesion_update.SetCellCellAdhesionEnergyParameter(0.5);
		adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.6);

     	TS_ASSERT_DELTA(adhesion_update.GetCellCellAdhesionEnergyParameter(), 0.5, 1e-12);
		TS_ASSERT_DELTA(adhesion_update.GetCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-12);

		adhesion_update.SetCellCellAdhesionEnergyParameter(0.1);
		adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.2);

        // Test get adhesion methods
		CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0u);
		CellPtr p_cell_1 = cell_population.GetCellUsingLocationIndex(1u);
        TS_ASSERT_DELTA(adhesion_update.GetCellCellAdhesionEnergy(p_cell_0,p_cell_1), 0.1, 1e-12);
        TS_ASSERT_DELTA(adhesion_update.GetCellBoundaryAdhesionEnergy(p_cell_0), 0.2, 1e-12);

		// Test EvaluateHamiltonianContribution()

        double gamma_cell_cell = adhesion_update.GetCellCellAdhesionEnergyParameter();
        double gamma_cell_boundary = adhesion_update.GetCellBoundaryAdhesionEnergyParameter();

		// Both points lie within cell 0
		TS_ASSERT_THROWS_THIS(adhesion_update.EvaluateHamiltonianContribution(0, 1, cell_population),
		                      "The current node and target node must not be in the same element.");

		// Both points lie within cell 1
		TS_ASSERT_THROWS_THIS(adhesion_update.EvaluateHamiltonianContribution(8, 9, cell_population),
							  "The current node and target node must not be in the same element.");

		// Both points lie within cell medium
		TS_ASSERT_THROWS_THIS(adhesion_update.EvaluateHamiltonianContribution(2, 3, cell_population),
		                       "At least one of the current node or target node must be in an element.");

		// Current site in cell 0; target site in cell medium
		double contribution = adhesion_update.EvaluateHamiltonianContribution(5, 6, cell_population);
		TS_ASSERT_DELTA(contribution, 2.0*gamma_cell_boundary, 1e-6);

		// Current site in cell medium; target site in cell 0
		contribution = adhesion_update.EvaluateHamiltonianContribution(6, 5, cell_population);
		TS_ASSERT_DELTA(contribution, 2.0*gamma_cell_boundary - gamma_cell_cell, 1e-6);

		// Current site in cell 0; target site in cell 1
		contribution = adhesion_update.EvaluateHamiltonianContribution(5, 9, cell_population);
		TS_ASSERT_DELTA(contribution, gamma_cell_cell, 1e-6);
    }

    void TestArchiveAdhesionPottsUpdateRule() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "AdhesionPottsUpdateRule.arch";

        {
            AdhesionPottsUpdateRule<2> update_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            update_rule.SetCellCellAdhesionEnergyParameter(0.5);
            update_rule.SetCellBoundaryAdhesionEnergyParameter(0.6);

            // Serialize via pointer to most abstract class possible
            AbstractPottsUpdateRule<2>* const p_update_rule = &update_rule;
            output_arch << p_update_rule;
        }

        {
            AbstractPottsUpdateRule<2>* p_update_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_update_rule;

            // Test the member data
            TS_ASSERT_DELTA((static_cast<AdhesionPottsUpdateRule<2>*>(p_update_rule))->GetCellCellAdhesionEnergyParameter(), 0.5, 1e-6);
            TS_ASSERT_DELTA((static_cast<AdhesionPottsUpdateRule<2>*>(p_update_rule))->GetCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-6);

            // Tidy up
            delete p_update_rule;
        }
    }

    void TestDifferentialAdhesionPottsUpdateRuleMethods() throw (Exception)
    {
        // Create a simple 2D PottsMesh with 4 elements
        PottsMeshGenerator<2> generator(5, 2, 2, 4, 2, 2, 1, 1, 1, true); // last bool makes elements start in bottom left
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Label cells 0 and 1
        boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);
        cells[0]->AddCellProperty(p_label);
        cells[1]->AddCellProperty(p_label);

        // Create cell population
        PottsBasedCellPopulation cell_population(*p_mesh, cells);

        // Create an update law system
        DifferentialAdhesionPottsUpdateRule<2> differential_adhesion_update;

        // Test get/set methods
        TS_ASSERT_DELTA(differential_adhesion_update.GetLabelledCellLabelledCellAdhesionEnergyParameter(), 0.1, 1e-12);
        TS_ASSERT_DELTA(differential_adhesion_update.GetLabelledCellCellAdhesionEnergyParameter(), 0.1, 1e-12);
        TS_ASSERT_DELTA(differential_adhesion_update.GetLabelledCellBoundaryAdhesionEnergyParameter(), 0.2, 1e-12);

        differential_adhesion_update.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.3);
        differential_adhesion_update.SetLabelledCellCellAdhesionEnergyParameter(0.4);
        differential_adhesion_update.SetLabelledCellBoundaryAdhesionEnergyParameter(0.5);

        TS_ASSERT_DELTA(differential_adhesion_update.GetLabelledCellLabelledCellAdhesionEnergyParameter(), 0.3, 1e-12);
        TS_ASSERT_DELTA(differential_adhesion_update.GetLabelledCellCellAdhesionEnergyParameter(), 0.4, 1e-12);
        TS_ASSERT_DELTA(differential_adhesion_update.GetLabelledCellBoundaryAdhesionEnergyParameter(), 0.5, 1e-12);

        // Test get adhesion methods
        CellPtr p_labelled_cell_0 = cell_population.GetCellUsingLocationIndex(0u);
        CellPtr p_labelled_cell_1 = cell_population.GetCellUsingLocationIndex(1u);
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(2u);
        CellPtr p_cell_1 = cell_population.GetCellUsingLocationIndex(3u);
        TS_ASSERT_DELTA(differential_adhesion_update.GetCellCellAdhesionEnergy(p_labelled_cell_0,p_labelled_cell_1), 0.3, 1e-12);
        TS_ASSERT_DELTA(differential_adhesion_update.GetCellCellAdhesionEnergy(p_labelled_cell_0,p_cell_0), 0.4, 1e-12);
        TS_ASSERT_DELTA(differential_adhesion_update.GetCellCellAdhesionEnergy(p_cell_0,p_cell_1), 0.1, 1e-12);// Default value from AdhesionPottsUpdateRule
        TS_ASSERT_DELTA(differential_adhesion_update.GetCellBoundaryAdhesionEnergy(p_labelled_cell_0), 0.5, 1e-12);
        TS_ASSERT_DELTA(differential_adhesion_update.GetCellBoundaryAdhesionEnergy(p_cell_0), 0.2, 1e-12); // Default value from AdhesionPottsUpdateRule

        // Set the parameters to facilitate calculations below
        differential_adhesion_update.SetCellCellAdhesionEnergyParameter(0.1);
        differential_adhesion_update.SetLabelledCellCellAdhesionEnergyParameter(0.2);
        differential_adhesion_update.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.3);
        differential_adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.4);
        differential_adhesion_update.SetLabelledCellBoundaryAdhesionEnergyParameter(0.5);

		// Test EvaluateHamiltonianContribution()

        double gamma_cell_0_cell_0 = differential_adhesion_update.GetCellCellAdhesionEnergyParameter();
        double gamma_cell_0_cell_1 = differential_adhesion_update.GetLabelledCellCellAdhesionEnergyParameter();
        double gamma_cell_1_cell_1 = differential_adhesion_update.GetLabelledCellLabelledCellAdhesionEnergyParameter();
        double gamma_cell_0_boundary = differential_adhesion_update.GetCellBoundaryAdhesionEnergyParameter();
        double gamma_cell_1_boundary = differential_adhesion_update.GetLabelledCellBoundaryAdhesionEnergyParameter();

		// Both points lie within cell 0
		TS_ASSERT_THROWS_THIS(differential_adhesion_update.EvaluateHamiltonianContribution(0, 1, cell_population),
		                      "The current node and target node must not be in the same element.");

		// Both points lie within cell medium
		TS_ASSERT_THROWS_THIS(differential_adhesion_update.EvaluateHamiltonianContribution(4, 9, cell_population),
		                       "At least one of the current node or target node must be in an element.");

		// Current site in cell 1; target site in cell medium
		double contribution = differential_adhesion_update.EvaluateHamiltonianContribution(8, 9, cell_population);
		TS_ASSERT_DELTA(contribution, gamma_cell_1_boundary, 1e-6);

		// Current site in cell 3; target site in cell medium
		contribution = differential_adhesion_update.EvaluateHamiltonianContribution(13, 14, cell_population);
		TS_ASSERT_DELTA(contribution, gamma_cell_0_boundary, 1e-6);

		// Current site in cell 0; target site in cell 1; both are labelled
		contribution = differential_adhesion_update.EvaluateHamiltonianContribution(1, 2, cell_population);
		TS_ASSERT_DELTA(contribution, gamma_cell_1_cell_1, 1e-6);

		// Current site in cell 2; target site in cell 3; both are unlabelled
		contribution = differential_adhesion_update.EvaluateHamiltonianContribution(11, 12, cell_population);
		TS_ASSERT_DELTA(contribution, gamma_cell_0_cell_0, 1e-6);

		// Current site in cell 0; target site in cell 2;
		contribution = differential_adhesion_update.EvaluateHamiltonianContribution(6, 11, cell_population);
		TS_ASSERT_DELTA(contribution, 2.0*gamma_cell_0_cell_1 - gamma_cell_0_cell_0, 1e-6);
    }

    void TestDifferentialArchiveAdhesionPottsUpdateRule() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "AdhesionPottsUpdateRule.arch";

        {
            DifferentialAdhesionPottsUpdateRule<2> update_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            update_rule.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.3);
            update_rule.SetLabelledCellCellAdhesionEnergyParameter(0.4);
            update_rule.SetCellCellAdhesionEnergyParameter(0.5);
            update_rule.SetLabelledCellBoundaryAdhesionEnergyParameter(0.6);
            update_rule.SetCellBoundaryAdhesionEnergyParameter(0.7);

            // Serialize via pointer to most abstract class possible
            AbstractPottsUpdateRule<2>* const p_update_rule = &update_rule;
            output_arch << p_update_rule;
        }

        {
            AbstractPottsUpdateRule<2>* p_update_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_update_rule;

            // Test the member data
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionPottsUpdateRule<2>*>(p_update_rule))->GetLabelledCellLabelledCellAdhesionEnergyParameter(), 0.3, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionPottsUpdateRule<2>*>(p_update_rule))->GetLabelledCellCellAdhesionEnergyParameter(), 0.4, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionPottsUpdateRule<2>*>(p_update_rule))->GetCellCellAdhesionEnergyParameter(), 0.5, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionPottsUpdateRule<2>*>(p_update_rule))->GetLabelledCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionPottsUpdateRule<2>*>(p_update_rule))->GetCellBoundaryAdhesionEnergyParameter(), 0.7, 1e-6);

            // Tidy up
            delete p_update_rule;
        }
    }


	void TestUpdateRuleOutputUpdateRuleInfo()
	{
		std::string output_directory = "TestPottsUpdateRulesOutputParameters";
		OutputFileHandler output_file_handler(output_directory, false);

		// Test with VolumeConstraintPottsUpdateRule
		VolumeConstraintPottsUpdateRule<2> volume_constraint;
        volume_constraint.SetDeformationEnergyParameter(0.1);
        volume_constraint.SetMatureCellTargetVolume(20);

		TS_ASSERT_EQUALS(volume_constraint.GetIdentifier(), "VolumeConstraintPottsUpdateRule-2");

		out_stream volume_constraint_parameter_file = output_file_handler.OpenOutputFile("volume_constraint_results.parameters");
		volume_constraint.OutputUpdateRuleInfo(volume_constraint_parameter_file);
		volume_constraint_parameter_file->close();

		std::string volume_constraint_results_dir = output_file_handler.GetOutputDirectoryFullPath();
		TS_ASSERT_EQUALS(system(("diff " + volume_constraint_results_dir + "volume_constraint_results.parameters notforrelease_cell_based/test/data/TestPottsUpdateRules/volume_constraint_results.parameters").c_str()), 0);

        // Test with SurfaceAreaConstraintPottsUpdateRule
        SurfaceAreaConstraintPottsUpdateRule<2> surface_area_constraint;
        surface_area_constraint.SetDeformationEnergyParameter(0.1);
        surface_area_constraint.SetMatureCellTargetSurfaceArea(20);

        TS_ASSERT_EQUALS(surface_area_constraint.GetIdentifier(), "SurfaceAreaConstraintPottsUpdateRule-2");

        out_stream surface_area_constraint_parameter_file = output_file_handler.OpenOutputFile("surface_area_constraint_results.parameters");
        surface_area_constraint.OutputUpdateRuleInfo(surface_area_constraint_parameter_file);
        surface_area_constraint_parameter_file->close();

        std::string surface_area_constraint_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + surface_area_constraint_results_dir + "surface_area_constraint_results.parameters notforrelease_cell_based/test/data/TestPottsUpdateRules/surface_area_constraint_results.parameters").c_str()), 0);

		// Test with AdhesionPottsUpdateRule
		AdhesionPottsUpdateRule<2> adhesion_update;
        adhesion_update.SetCellCellAdhesionEnergyParameter(0.3);
        adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.4);

		TS_ASSERT_EQUALS(adhesion_update.GetIdentifier(), "AdhesionPottsUpdateRule-2");

		out_stream adhesion_update_parameter_file = output_file_handler.OpenOutputFile("adhesion_update_results.parameters");
		adhesion_update.OutputUpdateRuleInfo(adhesion_update_parameter_file);
		adhesion_update_parameter_file->close();

		std::string adhesion_update_results_dir = output_file_handler.GetOutputDirectoryFullPath();
		TS_ASSERT_EQUALS(system(("diff " + adhesion_update_results_dir + "adhesion_update_results.parameters notforrelease_cell_based/test/data/TestPottsUpdateRules/adhesion_update_results.parameters").c_str()), 0);

        // Test with VolumeConstraintPottsUpdateRule
        DifferentialAdhesionPottsUpdateRule<2> differential_adhesion_update;
        differential_adhesion_update.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.3);
        differential_adhesion_update.SetLabelledCellCellAdhesionEnergyParameter(0.4);
        differential_adhesion_update.SetCellCellAdhesionEnergyParameter(0.5);
        differential_adhesion_update.SetLabelledCellBoundaryAdhesionEnergyParameter(0.6);
        differential_adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.7);

        TS_ASSERT_EQUALS(differential_adhesion_update.GetIdentifier(), "DifferentialAdhesionPottsUpdateRule-2");

        out_stream differential_adhesion_update_parameter_file = output_file_handler.OpenOutputFile("differential_adhesion_update_results.parameters");
        differential_adhesion_update.OutputUpdateRuleInfo(differential_adhesion_update_parameter_file);
        differential_adhesion_update_parameter_file->close();

        std::string differential_adhesion_update_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + differential_adhesion_update_results_dir + "differential_adhesion_update_results.parameters notforrelease_cell_based/test/data/TestPottsUpdateRules/differential_adhesion_update_results.parameters").c_str()), 0);
	}
};

#endif /*TESTPOTTSUPDATERULES_HPP_*/
