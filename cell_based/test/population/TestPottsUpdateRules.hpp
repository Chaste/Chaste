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
#include "ChemotaxisPottsUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "PottsMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestPottsUpdateRules : public AbstractCellBasedTestSuite
{
public:

    void TestVolumeConstraintPottsUpdateRuleIn2d()
    {
        // Create a simple 2D PottsMesh with 2 elements
        PottsMeshGenerator<2> generator(4, 1, 2, 4, 2, 2, 1, 1, 1, true); // last bool makes elements start in bottom left
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

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

    void TestVolumeConstraintPottsUpdateRuleIn3d()
    {
        // Create a simple 2D PottsMesh with 2 elements
        PottsMeshGenerator<3> generator(4, 2, 2, 2, 1, 2, 4, 1, 2, true);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 32u);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);

        // Create an update law system
        VolumeConstraintPottsUpdateRule<3> volume_constraint;

        // Test EvaluateHamiltonianContribution()
        volume_constraint.SetDeformationEnergyParameter(0.5);
        volume_constraint.SetMatureCellTargetVolume(8.0);

        double alpha = volume_constraint.GetDeformationEnergyParameter();

        // Both points lie within cell 0
        TS_ASSERT_THROWS_THIS(volume_constraint.EvaluateHamiltonianContribution(0, 1, cell_population),
                             "The current node and target node must not be in the same element.");

        // Both points lie within cell 1
        TS_ASSERT_THROWS_THIS(volume_constraint.EvaluateHamiltonianContribution(2, 3, cell_population),
                             "The current node and target node must not be in the same element.");

        // Both points lie within cell medium
        TS_ASSERT_THROWS_THIS(volume_constraint.EvaluateHamiltonianContribution(16, 17, cell_population),
                             "At least one of the current node or target node must be in an element.");

        // Current site in cell 0; target site in cell medium
        double contribution = volume_constraint.EvaluateHamiltonianContribution(9, 17, cell_population);
        TS_ASSERT_DELTA(contribution, alpha, 1e-6);

        // Current site in cell medium; target site in cell 0
        contribution = volume_constraint.EvaluateHamiltonianContribution(17, 9, cell_population);
        TS_ASSERT_DELTA(contribution, alpha, 1e-6);

        // Current site in cell 0; target site in cell 1
        contribution = volume_constraint.EvaluateHamiltonianContribution(9, 10, cell_population);
        TS_ASSERT_DELTA(contribution, 2.0*alpha, 1e-6);
    }

    void TestArchiveVolumeConstraintPottsUpdateRule()
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

    void TestSurfaceAreaConstraintPottsUpdateRuleIn2d()
    {
        // Create a simple 2D PottsMesh with 2 elements
        PottsMeshGenerator<2> generator(4, 1, 2, 4, 2, 2, 1, 1, 1, true); // last bool makes elements start in bottom left
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

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
        TS_ASSERT_DELTA(contribution, 2.0*2.0*gamma, 1e-6);

        // Current site in cell 0; target site in cell medium
        contribution = surface_area_constraint.EvaluateHamiltonianContribution(5, 6, cell_population);
        TS_ASSERT_DELTA(contribution, 2.0*2.0*gamma, 1e-6);

        // Current site in cell medium; target site in cell 0
        contribution = surface_area_constraint.EvaluateHamiltonianContribution(6, 5, cell_population);
        TS_ASSERT_DELTA(contribution, 0.0, 1e-6);

        // Current site in cell 0; target site in cell 1
        contribution = surface_area_constraint.EvaluateHamiltonianContribution(5, 9, cell_population);
        TS_ASSERT_DELTA(contribution, 2.0*2.0*gamma, 1e-6);

        // Current site in cell 0; target site in cell 1 (diagonal neighbour)
        contribution = surface_area_constraint.EvaluateHamiltonianContribution(5, 8, cell_population);
        TS_ASSERT_DELTA(contribution, 2.0*2.0*gamma, 1e-6);
    }

    void TestSurfaceAreaConstraintPottsUpdateRuleIn3d()
    {
        // Create a simple 3D PottsMesh with 2 elements
        PottsMeshGenerator<3> generator(4, 2, 2, 4, 1, 2, 4, 1, 2, true);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 64u);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);

        // Create an update law system
        SurfaceAreaConstraintPottsUpdateRule<3> surface_area_constraint;

        surface_area_constraint.SetDeformationEnergyParameter(0.5);
        surface_area_constraint.SetMatureCellTargetSurfaceArea(24.0);

        // Test EvaluateHamiltonianContribution()
        double gamma = surface_area_constraint.GetDeformationEnergyParameter();

        // Both points lie within cell 0
        TS_ASSERT_THROWS_THIS(surface_area_constraint.EvaluateHamiltonianContribution(0, 1, cell_population),
                              "The current node and target node must not be in the same element.");

        // Both points lie within cell 1
        TS_ASSERT_THROWS_THIS(surface_area_constraint.EvaluateHamiltonianContribution(2, 3, cell_population),
                              "The current node and target node must not be in the same element.");

        // Both points lie within cell medium
        TS_ASSERT_THROWS_THIS(surface_area_constraint.EvaluateHamiltonianContribution(8, 9, cell_population),
                               "At least one of the current node or target node must be in an element.");

        // Current site in cell 0; target site in cell medium Cell on edge of domain
        double contribution = surface_area_constraint.EvaluateHamiltonianContribution(4, 8, cell_population);
        TS_ASSERT_DELTA(contribution, 4.0*4.0*gamma, 1e-6);

        // Current site in cell medium; target site in cell 1
        contribution = surface_area_constraint.EvaluateHamiltonianContribution(10, 6, cell_population);
        TS_ASSERT_DELTA(contribution, 0.0, 1e-6);

        // Current site in cell 0; target site in cell 1
        contribution = surface_area_constraint.EvaluateHamiltonianContribution(17, 18, cell_population);
        TS_ASSERT_DELTA(contribution, 4.0*4.0*gamma, 1e-6);

        // Current site in cell 0; target site in cell 1 (diagonal switch)
        contribution = surface_area_constraint.EvaluateHamiltonianContribution(1, 18, cell_population);
        TS_ASSERT_DELTA(contribution, 4.0*4.0*gamma, 1e-6);

        // Current site in cell 0; target site in medium
        contribution = surface_area_constraint.EvaluateHamiltonianContribution(5, 9, cell_population);
        TS_ASSERT_DELTA(contribution, 4.0*4.0*gamma, 1e-6);
    }

    void TestArchiveSurfaceAreaConstraintPottsUpdateRule()
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

    void TestAdhesionPottsUpdateRuleIn2d()
    {
        // Create a simple 2D PottsMesh with 2 elements
        PottsMeshGenerator<2> generator(4, 1, 2, 4, 2, 2, 1, 1, 1, true); // last bool makes elements start in bottom left
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

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

        // Test EvaluateHamiltonianContribution(). Note that the boundary of the domain is a void not medium.

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

    void TestAdhesionPottsUpdateRuleIn3d()
    {
        // Create a simple 2D PottsMesh with 2 elements
        PottsMeshGenerator<3> generator(4, 2, 2, 2, 1, 2, 4, 1, 2, true);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 32u);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Label cells 0 and 1
        MAKE_PTR(CellLabel, p_label);
        cells[0]->AddCellProperty(p_label);
        cells[1]->AddCellProperty(p_label);

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);

        // Create an update law system
        AdhesionPottsUpdateRule<3> adhesion_update;

        // Set the parameters to facilitate calculations below.
        adhesion_update.SetCellCellAdhesionEnergyParameter(0.1);
        adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.2);

        // Test EvaluateHamiltonianContribution(). Note that the boundary of the domain is a void not medium.

        double gamma_cell_cell = adhesion_update.GetCellCellAdhesionEnergyParameter();
        double gamma_cell_boundary = adhesion_update.GetCellBoundaryAdhesionEnergyParameter();

        // Both points lie within cell 0
        TS_ASSERT_THROWS_THIS(adhesion_update.EvaluateHamiltonianContribution(0, 1, cell_population),
                              "The current node and target node must not be in the same element.");

        // Both points lie within cell 1
        TS_ASSERT_THROWS_THIS(adhesion_update.EvaluateHamiltonianContribution(2, 3, cell_population),
                              "The current node and target node must not be in the same element.");

        // Both points lie within cell medium
        TS_ASSERT_THROWS_THIS(adhesion_update.EvaluateHamiltonianContribution(16, 17, cell_population),
                               "At least one of the current node or target node must be in an element.");

        // Current site in cell 0; target site in cell medium
        double contribution = adhesion_update.EvaluateHamiltonianContribution(9, 17, cell_population);
        TS_ASSERT_DELTA(contribution, 3.0*gamma_cell_boundary, 1e-6);

        // Current site in cell medium; target site in cell 0
        contribution = adhesion_update.EvaluateHamiltonianContribution(17, 9, cell_population);
        TS_ASSERT_DELTA(contribution, 3.0*gamma_cell_boundary - gamma_cell_cell, 1e-6);

        // Current site in cell 0; target site in cell 1
        contribution = adhesion_update.EvaluateHamiltonianContribution(9, 10, cell_population);
        TS_ASSERT_DELTA(contribution, 2.0*gamma_cell_cell, 1e-6);
    }

    void TestArchiveAdhesionPottsUpdateRule()
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

    void TestDifferentialAdhesionPottsUpdateRuleIn2d()
    {
        // Create a simple 2D PottsMesh with 4 elements
        PottsMeshGenerator<2> generator(5, 2, 2, 4, 2, 2, 1, 1, 1, true); // last bool makes elements start in bottom left
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Label cells 0 and 1
        MAKE_PTR(CellLabel, p_label);
        cells[0]->AddCellProperty(p_label);
        cells[1]->AddCellProperty(p_label);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

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

        // Test EvaluateHamiltonianContribution(). Note that the boundary of the domain is a void not medium

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

    void TestDifferentialAdhesionPottsUpdateRuleIn3d()
    {
        // Create a simple 2D PottsMesh with 4 elements
        PottsMeshGenerator<3> generator(4, 2, 2, 2, 2, 1, 4, 1, 2, true);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 32u);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Label cells 0 and 1
        MAKE_PTR(CellLabel, p_label);
        cells[0]->AddCellProperty(p_label);
        cells[1]->AddCellProperty(p_label);

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);

        // Create an update law system
        DifferentialAdhesionPottsUpdateRule<3> differential_adhesion_update;

        // Set the parameters to facilitate calculations below
        differential_adhesion_update.SetCellCellAdhesionEnergyParameter(0.1);
        differential_adhesion_update.SetLabelledCellCellAdhesionEnergyParameter(0.2);
        differential_adhesion_update.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.3);
        differential_adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.4);
        differential_adhesion_update.SetLabelledCellBoundaryAdhesionEnergyParameter(0.5);

        // Test EvaluateHamiltonianContribution(). Note that the boundary of the domain is a void not medium.

        double gamma_cell_0_cell_0 = differential_adhesion_update.GetCellCellAdhesionEnergyParameter();
        double gamma_cell_0_cell_1 = differential_adhesion_update.GetLabelledCellCellAdhesionEnergyParameter();
        double gamma_cell_1_cell_1 = differential_adhesion_update.GetLabelledCellLabelledCellAdhesionEnergyParameter();
        double gamma_cell_0_boundary = differential_adhesion_update.GetCellBoundaryAdhesionEnergyParameter();
        double gamma_cell_1_boundary = differential_adhesion_update.GetLabelledCellBoundaryAdhesionEnergyParameter();

        // Both points lie within cell 0
        TS_ASSERT_THROWS_THIS(differential_adhesion_update.EvaluateHamiltonianContribution(0, 1, cell_population),
                              "The current node and target node must not be in the same element.");

        // Both points lie within cell medium
        TS_ASSERT_THROWS_THIS(differential_adhesion_update.EvaluateHamiltonianContribution(20, 21, cell_population),
                               "At least one of the current node or target node must be in an element.");

        // Current site in cell 1; target site in cell medium
        double contribution = differential_adhesion_update.EvaluateHamiltonianContribution(9, 17, cell_population);
        TS_ASSERT_DELTA(contribution, 3.0*gamma_cell_1_boundary, 1e-6);

        // Current site in cell 3; target site in cell medium
        contribution = differential_adhesion_update.EvaluateHamiltonianContribution(12, 20, cell_population);
        TS_ASSERT_DELTA(contribution, 2.0 * gamma_cell_0_boundary, 1e-6);

        // Current site in cell 0; target site in cell 1; both are labelled
        contribution = differential_adhesion_update.EvaluateHamiltonianContribution(9, 10, cell_population);
        TS_ASSERT_DELTA(contribution, gamma_cell_1_cell_1, 1e-6);

        // Current site in cell 2; target site in cell 3; both are unlabelled
        contribution = differential_adhesion_update.EvaluateHamiltonianContribution(13, 14, cell_population);
        TS_ASSERT_DELTA(contribution, gamma_cell_0_cell_0, 1e-6);

        // Current site in cell 0; target site in cell 2;
        contribution = differential_adhesion_update.EvaluateHamiltonianContribution(9, 13, cell_population);
        TS_ASSERT_DELTA(contribution, 2.0*gamma_cell_0_cell_1 - gamma_cell_0_cell_0 + gamma_cell_1_boundary - gamma_cell_0_boundary, 1e-6);
    }

    void TestArchiveDifferentialAdhesionPottsUpdateRule()
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

    void TestChemotaxisPottsUpdateRuleIn2d()
    {
        // Create a simple 2D PottsMesh with 2 elements
        PottsMeshGenerator<2> generator(4, 1, 2, 4, 1, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an update law system
        ChemotaxisPottsUpdateRule<2> chemotaxis_update;

        // Test EvaluateHamiltonianContribution()

        // target site above current site
        double contribution = chemotaxis_update.EvaluateHamiltonianContribution(10, 14, cell_population);
        TS_ASSERT_DELTA(contribution, -0.2, 1e-6);

        // target site below current site
        contribution = chemotaxis_update.EvaluateHamiltonianContribution(6, 2, cell_population);
        TS_ASSERT_DELTA(contribution, 0.2, 1e-6);

        // target site diagonally above current site (to right)
        contribution = chemotaxis_update.EvaluateHamiltonianContribution(10, 15, cell_population);
        TS_ASSERT_DELTA(contribution, -0.4, 1e-6);

        // target site diagonally above current site (to left)
        contribution = chemotaxis_update.EvaluateHamiltonianContribution(10, 13, cell_population);
        TS_ASSERT_DELTA(contribution, 0.0, 1e-6);
    }

    void TestArchiveChemotaxisPottsUpdateRule()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ChemotaxisPottsUpdateRule.arch";

        {
            ChemotaxisPottsUpdateRule<2> update_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            // Currently none to test

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
            // Currently none to test

            // Tidy up
            delete p_update_rule;
        }
    }

    void TestUpdateRuleOutputUpdateRuleInfo()
    {
        EXIT_IF_PARALLEL;
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
        FileComparison( volume_constraint_results_dir + "volume_constraint_results.parameters", "cell_based/test/data/TestPottsUpdateRules/volume_constraint_results.parameters").CompareFiles();

        // Test with SurfaceAreaConstraintPottsUpdateRule
        SurfaceAreaConstraintPottsUpdateRule<2> surface_area_constraint;
        surface_area_constraint.SetDeformationEnergyParameter(0.1);
        surface_area_constraint.SetMatureCellTargetSurfaceArea(20);

        TS_ASSERT_EQUALS(surface_area_constraint.GetIdentifier(), "SurfaceAreaConstraintPottsUpdateRule-2");

        out_stream surface_area_constraint_parameter_file = output_file_handler.OpenOutputFile("surface_area_constraint_results.parameters");
        surface_area_constraint.OutputUpdateRuleInfo(surface_area_constraint_parameter_file);
        surface_area_constraint_parameter_file->close();

        std::string surface_area_constraint_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        FileComparison( surface_area_constraint_results_dir + "surface_area_constraint_results.parameters", "cell_based/test/data/TestPottsUpdateRules/surface_area_constraint_results.parameters").CompareFiles();

        // Test with AdhesionPottsUpdateRule
        AdhesionPottsUpdateRule<2> adhesion_update;
        adhesion_update.SetCellCellAdhesionEnergyParameter(0.3);
        adhesion_update.SetCellBoundaryAdhesionEnergyParameter(0.4);

        TS_ASSERT_EQUALS(adhesion_update.GetIdentifier(), "AdhesionPottsUpdateRule-2");

        out_stream adhesion_update_parameter_file = output_file_handler.OpenOutputFile("adhesion_update_results.parameters");
        adhesion_update.OutputUpdateRuleInfo(adhesion_update_parameter_file);
        adhesion_update_parameter_file->close();

        std::string adhesion_update_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        FileComparison( adhesion_update_results_dir + "adhesion_update_results.parameters", "cell_based/test/data/TestPottsUpdateRules/adhesion_update_results.parameters").CompareFiles();

        // Test with DifferentialAdhesionPottsUpdateRule
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
        FileComparison(differential_adhesion_update_results_dir + "differential_adhesion_update_results.parameters",
                       "cell_based/test/data/TestPottsUpdateRules/differential_adhesion_update_results.parameters").CompareFiles();

        // Test with DifferentialAdhesionPottsUpdateRule
        ChemotaxisPottsUpdateRule<2> chemotaxis_update;

        TS_ASSERT_EQUALS(chemotaxis_update.GetIdentifier(), "ChemotaxisPottsUpdateRule-2");

        out_stream chemotaxis_update_parameter_file = output_file_handler.OpenOutputFile("chemotaxis_update_results.parameters");
        chemotaxis_update.OutputUpdateRuleInfo(chemotaxis_update_parameter_file);
        chemotaxis_update_parameter_file->close();

        std::string chemotaxis_update_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        FileComparison(chemotaxis_update_results_dir + "chemotaxis_update_results.parameters",
                       "cell_based/test/data/TestPottsUpdateRules/chemotaxis_update_results.parameters").CompareFiles();
    }
};

#endif /*TESTPOTTSUPDATERULES_HPP_*/
