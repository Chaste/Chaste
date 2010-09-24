/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef TESTFORCES2NOTFORRELEASE_HPP_
#define TESTFORCES2NOTFORRELEASE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WntConcentration.hpp"
#include "CryptProjectionForce.hpp"
#include "HoneycombMutableVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
//#include "LinearSpringWithVariableSpringConstantsForce.hpp"
//#include "CellwiseDataGradient.hpp"
//#include "VertexCryptBoundaryForce.hpp"
//#include "ApcTwoHitCellMutationState.hpp"
//#include "WildTypeCellMutationState.hpp"
//#include "CellLabel.hpp"

class TestForces2NotForRelease : public AbstractCellBasedTestSuite
{
public:
    // This test contains cells of 2 mutation types, wildtype and labelled type,
    // on a larger mesh so that we can test interaction of 2 cells with labelled types.
    // It asserts that neighboring cells have the correct adhesion parameter for difference
    // pairs of nodes.
    void TestNagaiHondaForceForCellsWithTwoMutationTypes() throw (Exception)
    {
        // Create a simple 2D MutableVertexMesh with four cells
        HoneycombMutableVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMutableMesh();

        // Set up cells.
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);

        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(TRANSIT);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -2.0;
            p_cell->SetBirthTime(birth_time);

            if (elem_index == 0 || elem_index == 2) // cells chosen for coverage
            {
                p_cell->AddCellProperty(p_label);
            }
            cells.push_back(p_cell);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force system
        NagaiHondaDifferentialAdhesionForce<2> force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&force);

        // Check mutation state
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_EQUALS(cells[i]->GetMutationState()->IsType<WildTypeCellMutationState>(), true);
            if (i==0 || i==2)
            {
                TS_ASSERT_EQUALS(cells[i]->HasCellProperty<CellLabel>(), true);
            }
        }

        // There are two combinations of type WILD_WILD
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(7), p_mesh->GetNode(9), cell_population), WILD_WILD);
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(9), p_mesh->GetNode(7), cell_population), WILD_WILD);
        TS_ASSERT_DELTA(force.GetAdhesionParameterDifferentialAddition(p_mesh->GetNode(9), p_mesh->GetNode(7), WILD_WILD), 0.01, 1e-4);

        // There are two combinations of type WILD_LABELLED
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(6), p_mesh->GetNode(9), cell_population), WILD_LABELLED);
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(9), p_mesh->GetNode(6), cell_population), WILD_LABELLED);
        TS_ASSERT_DELTA(force.GetAdhesionParameterDifferentialAddition(p_mesh->GetNode(9), p_mesh->GetNode(6), WILD_LABELLED), 1.0, 1e-4);

        // There are two combinations of type LABELLED_LABELLED
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(6), p_mesh->GetNode(8), cell_population), LABELLED_LABELLED);
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(8), p_mesh->GetNode(6), cell_population), LABELLED_LABELLED);
        TS_ASSERT_DELTA(force.GetAdhesionParameterDifferentialAddition(p_mesh->GetNode(9), p_mesh->GetNode(7), LABELLED_LABELLED), 0.01, 1e-4);

        // There is one combination of type OTHER (labelled/void)
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(0), p_mesh->GetNode(3), cell_population), OTHER);
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(3), p_mesh->GetNode(0), cell_population), OTHER);

        // There is one combination of type OTHER (wild type/void)
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(10), p_mesh->GetNode(13), cell_population), OTHER);
        TS_ASSERT_EQUALS(force.GetCombinationCellTypes(p_mesh->GetNode(13), p_mesh->GetNode(10), cell_population), OTHER);

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        force.AddForceContribution(node_forces, cell_population);

        // Check some example forces (these will change if you modify the adhesion parameters)
        TS_ASSERT_DELTA(node_forces[0][0], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_forces[0][1], -14.0135, 1e-4);

        TS_ASSERT_DELTA(node_forces[10][0], 12.1361, 1e-4);
        TS_ASSERT_DELTA(node_forces[10][1], -7.0067, 1e-4);
    }


    void TestArchivingNagaiHondaDifferentialAdhesionForce() throw (Exception)
    {
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "nagai_honda_differential_adhesion.arch";

        {
            // Construct a 2D vertex mesh consisting of a single element
            std::vector<Node<2>*> nodes;
            unsigned num_nodes = 20;
            std::vector<double> angles = std::vector<double>(num_nodes);

            for (unsigned i=0; i<num_nodes; i++)
            {
                angles[i] = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
                nodes.push_back(new Node<2>(i, true, cos(angles[i]), sin(angles[i])));
            }

            std::vector<VertexElement<2,2>*> elements;
            elements.push_back(new VertexElement<2,2>(0, nodes));

            double cell_swap_threshold = 0.01;
            double edge_division_threshold = 2.0;
            MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

            // Set up the cell
            std::vector<CellPtr> cells;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(-1.0);
            cells.push_back(p_cell);

            // Create cell population
            VertexBasedCellPopulation<2> cell_population(mesh, cells);
            cell_population.InitialiseCells();

            // Create a force system
            NagaiHondaDifferentialAdhesionForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            NagaiHondaDifferentialAdhesionForce<2>* const p_force = &force;
            output_arch << p_force;

            // Tidy up
            CellBasedConfig::Instance()->Reset();
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            NagaiHondaDifferentialAdhesionForce<2>* p_force;

            // Restore from the archive
            input_arch >> p_force;

            // Tidy up
            delete p_force;
        }
    }

    void TestForceOutputParameters()
    {
        std::string output_directory = "TestNotForReleaseForceOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with CryptProjectionForce
        CryptProjectionForce projection_force;
        TS_ASSERT_EQUALS(projection_force.GetIdentifier(), "CryptProjectionForce");

        out_stream projection_force_parameter_file = output_file_handler.OpenOutputFile("projection_results.parameters");
        projection_force.OutputForceParameters(projection_force_parameter_file);
        projection_force_parameter_file->close();

        std::string variable_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + variable_force_results_dir + "projection_results.parameters notforrelease_cell_based/test/data/TestForcesNotForRelease/projection_results.parameters").c_str()), 0);

        // Test with NagaiHondaForceDifferentialAdhesionForce
        NagaiHondaDifferentialAdhesionForce<2> differential_force;
        TS_ASSERT_EQUALS(differential_force.GetIdentifier(), "NagaiHondaDifferentialAdhesionForce-2");

        out_stream differential_force_parameter_file = output_file_handler.OpenOutputFile("differential_results.parameters");
        differential_force.OutputForceParameters(differential_force_parameter_file);
        differential_force_parameter_file->close();

        std::string differential_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + differential_force_results_dir + "differential_results.parameters notforrelease_cell_based/test/data/TestForcesNotForRelease/differential_results.parameters").c_str()), 0);
    }
};

#endif /*TESTFORCES2NOTFORRELEASE_HPP_*/

