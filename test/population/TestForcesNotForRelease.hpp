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
#ifndef TESTFORCESNOTFORRELEASE_HPP_
#define TESTFORCESNOTFORRELEASE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombMutableVertexMeshGenerator.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "CryptProjectionForce.hpp"
#include "CellwiseDataGradient.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "VertexCryptBoundaryForce.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WntConcentration.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"

class TestForcesNotForRelease : public AbstractCellBasedTestSuite
{
public:

    void TestCryptProjectionForceMethods() throw (Exception)
    {
        // Create a mesh
        unsigned num_cells_width = 10;
        unsigned num_cells_depth = 10;
        unsigned thickness_of_ghost_layer = 0;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Centre the mesh at (0,0)
        ChasteCuboid<2> bounding_box=p_mesh->CalculateBoundingBox();
        double width_of_mesh = (num_cells_width/(num_cells_width+2.0*thickness_of_ghost_layer))*(bounding_box.GetWidth(0));
        double height_of_mesh = (num_cells_depth/(num_cells_depth+2.0*thickness_of_ghost_layer))*(bounding_box.GetWidth(1));

        p_mesh->Translate(-width_of_mesh/2, -height_of_mesh/2);

        // Create some cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            CellPtr p_cell(new Cell(p_state, p_model));

            if (i==4 || i==5)
            {
                p_cell->SetBirthTime(-0.5);
            }
            else
            {
                p_cell->SetBirthTime(-10.0);
            }
            cells.push_back(p_cell);
        }

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        std::set<CellPtr> cell_pair_4_5 = cell_population.CreateCellPair(cell_population.GetCellUsingLocationIndex(4), cell_population.GetCellUsingLocationIndex(5));
        cell_population.MarkSpring(cell_pair_4_5);

        // Create a spring system with crypt surface z = 2*r
        WntConcentration<2>::Instance()->SetCryptProjectionParameterA(2.0);
        WntConcentration<2>::Instance()->SetCryptProjectionParameterB(1.0);
        CryptProjectionForce crypt_projection_force;

        // Test get methods
        TS_ASSERT_DELTA(crypt_projection_force.GetA(), 2.0, 1e-12);
        TS_ASSERT_DELTA(crypt_projection_force.GetB(), 1.0, 1e-12);

        // Test crypt height and gradient calculations
        c_vector<double, 2> node_location_2d = p_mesh->GetNode(0)->rGetLocation();
        TS_ASSERT_DELTA(crypt_projection_force.CalculateCryptSurfaceHeightAtPoint(node_location_2d), 2.0*pow(norm_2(node_location_2d),1.0), 1e-12);
        TS_ASSERT_DELTA(crypt_projection_force.CalculateCryptSurfaceDerivativeAtPoint(node_location_2d), 2.0, 1e-12);

        // Test updating of mNode3dLocationMap
        crypt_projection_force.UpdateNode3dLocationMap(cell_population);

        // Move a node slightly
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = node_location_2d[0]+0.05;
        new_point.rGetLocation()[1] = node_location_2d[1];
        p_mesh->SetNode(0, new_point, false);

        // Test UpdateNode3dLocationMap()

        c_vector<double, 2> new_node_location_2d;
        new_node_location_2d[0] = new_point.rGetLocation()[0];
        new_node_location_2d[1] = new_point.rGetLocation()[1];

        crypt_projection_force.UpdateNode3dLocationMap(cell_population);

        // Check the map updates correctly (note that we have used no ghost nodes, so the map does contain 0)
        c_vector<double, 3> calculated_new_node_location_3d = crypt_projection_force.mNode3dLocationMap[0];
        c_vector<double, 3> correct_new_node_location_3d;

        correct_new_node_location_3d[0] = new_node_location_2d[0];
        correct_new_node_location_3d[1] = new_node_location_2d[1];
        correct_new_node_location_3d[2] = crypt_projection_force.CalculateCryptSurfaceHeightAtPoint(new_node_location_2d);

        TS_ASSERT_DELTA(calculated_new_node_location_3d[0], correct_new_node_location_3d[0], 1e-12);
        TS_ASSERT_DELTA(calculated_new_node_location_3d[1], correct_new_node_location_3d[1], 1e-12);
        TS_ASSERT_DELTA(calculated_new_node_location_3d[2], correct_new_node_location_3d[2], 1e-12);

        // Test force calculation on a normal spring

        c_vector<double,2> force_on_spring; // between nodes 0 and 1

        // Find one of the elements that nodes 0 and 1 live on
        ChastePoint<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0] + 0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01;

        unsigned elem_index = p_mesh->GetContainingElementIndex(new_point2, false);
        Element<2,2>* p_element = p_mesh->GetElement(elem_index);

        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                            p_element->GetNodeGlobalIndex(0),
                                                                            cell_population);

        TS_ASSERT_DELTA(force_on_spring[0], -5.7594, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1],  0.0230, 1e-4);


        // Test force calculation with a cutoff

        double dist = norm_2(p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(),
                             p_element->GetNode(1)->rGetLocation()));

        crypt_projection_force.UseCutoffPoint(dist - 0.1);

        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                            p_element->GetNodeGlobalIndex(0),
                                                                            cell_population);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);

        // Test force calculation for a pair of newly born neighbouring cells
        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(4, 5, cell_population);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);

        cell_population.UnmarkSpring(cell_pair_4_5);

        // For coverage, test force calculation for a pair of neighbouring apoptotic cells
        cell_population.GetCellUsingLocationIndex(6)->StartApoptosis();
        cell_population.GetCellUsingLocationIndex(7)->StartApoptosis();
        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(6, 7, cell_population);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);

        // Test force calculation for a particular node

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }

        crypt_projection_force.AddForceContribution(node_forces, cell_population);

        TS_ASSERT_DELTA(node_forces[0][0], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_forces[0][1], 0.0, 1e-4);

        // Test that in the case of a flat crypt surface (mA=mB=0), the results are the same as for Meineke2001SpringSystem
        WntConcentration<2>::Instance()->SetCryptProjectionParameterA(0.001);
        WntConcentration<2>::Instance()->SetCryptProjectionParameterB(0.001);
        CryptProjectionForce flat_crypt_projection_force;
        GeneralisedLinearSpringForce<2> linear_force;

        // Normally this would be set up at the start of rCalculateforcesOfEachNode
        flat_crypt_projection_force.UpdateNode3dLocationMap(cell_population);

        for (MeshBasedCellPopulation<2>::SpringIterator spring_iterator = cell_population.SpringsBegin();
            spring_iterator != cell_population.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

            c_vector<double, 2> force_flat = flat_crypt_projection_force.CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, cell_population);
            c_vector<double, 2> force_meineke = linear_force.CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, cell_population);

            TS_ASSERT_DELTA( force_flat[0], force_meineke[0], 1e-3);
            TS_ASSERT_DELTA( force_flat[1], force_meineke[1], 1e-3);
        }

        WntConcentration<2>::Destroy();
    }

    /**
     * Note: WntBasedChemotaxis should be possible in other force laws. If/when
     * this is implemented, this test should be moved to somewhere more appropriate.
     */
    void TestCryptProjectionForceWithWntBasedChemotaxis() throw (Exception)
    {
        // Create a mesh
        unsigned num_cells_width = 10;
        unsigned num_cells_depth = 10;
        unsigned thickness_of_ghost_layer = 0;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Centre the mesh at (0,0)
        ChasteCuboid<2> bounding_box=p_mesh->CalculateBoundingBox();
        double width_of_mesh = (num_cells_width/(num_cells_width+2.0*thickness_of_ghost_layer))*(bounding_box.GetWidth(0));
        double height_of_mesh = (num_cells_depth/(num_cells_depth+2.0*thickness_of_ghost_layer))*(bounding_box.GetWidth(1));

        p_mesh->Translate(-width_of_mesh/2, -height_of_mesh/2);

        // Create some cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(-10.0);
            cells.push_back(p_cell);
        }

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        std::set<CellPtr> cell_pair_4_5 = cell_population.CreateCellPair(cell_population.GetCellUsingLocationIndex(4), cell_population.GetCellUsingLocationIndex(5));
        cell_population.MarkSpring(cell_pair_4_5);

        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);

        // Create a spring system with crypt surface z = 2*r
        WntConcentration<2>::Instance()->SetCryptProjectionParameterA(2.0);
        WntConcentration<2>::Instance()->SetCryptProjectionParameterB(1.0);
        CryptProjectionForce crypt_projection_force;

        crypt_projection_force.SetWntChemotaxis(false);

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > old_node_forces;
        old_node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             old_node_forces.push_back(zero_vector<double>(2));
        }

        // Calculate node forces
        crypt_projection_force.AddForceContribution(old_node_forces, cell_population);

        // Store the force of a particular node without Wnt-chemotaxis
        c_vector<double,2> old_force = old_node_forces[11];

        // Now turn on Wnt-chemotaxis
        crypt_projection_force.SetWntChemotaxis(true);

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > new_node_forces;
        new_node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             new_node_forces.push_back(zero_vector<double>(2));
        }

        // Calculate node forces
        crypt_projection_force.AddForceContribution(new_node_forces, cell_population);

        // Store the force of the same node, but now with Wnt-chemotaxis
        c_vector<double,2> new_force = new_node_forces[11];

        double wnt_chemotaxis_strength = crypt_projection_force.GetWntChemotaxisStrength();
        c_vector<double,2> wnt_component = wnt_chemotaxis_strength*WntConcentration<2>::Instance()->GetWntGradient(cells[11]);

        TS_ASSERT_DELTA(new_force[0], old_force[0]+wnt_component[0], 1e-4);
        TS_ASSERT_DELTA(new_force[1], old_force[1]+wnt_component[1], 1e-4);

        WntConcentration<2>::Destroy();
    }

    void TestCryptProjectionForceWithArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "crypt_projection_spring_system.arch";

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");

            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

            std::vector<CellPtr> cells;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
                p_model->SetCellProliferativeType(STEM);

                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetBirthTime(-50.0);
                cells.push_back(p_cell);
            }

            MeshBasedCellPopulation<2> crypt(mesh, cells);
            WntConcentration<2>::Instance()->SetCryptProjectionParameterA(1.0);
            WntConcentration<2>::Instance()->SetCryptProjectionParameterB(2.0);

            // Create force object
            CryptProjectionForce crypt_projection_force;

            TS_ASSERT_DELTA(crypt_projection_force.GetWntChemotaxisStrength(), 100.0, 1e-6);
            crypt_projection_force.SetWntChemotaxisStrength(15.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            CryptProjectionForce* const p_crypt_projection_force = &crypt_projection_force;

            p_crypt_projection_force->UseCutoffPoint(1.1);

            output_arch << p_crypt_projection_force;
            WntConcentration<2>::Destroy();
        }

        {
            ArchiveLocationInfo::SetMeshPathname("mesh/test/data/", "square_2_elements");

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            CryptProjectionForce* p_crypt_projection_force;

            // Restore from the archive
            input_arch >> p_crypt_projection_force;

            // Test the member data
            TS_ASSERT_EQUALS(p_crypt_projection_force->mUseCutoffPoint, true);
            TS_ASSERT_DELTA(p_crypt_projection_force->GetA(), 1.0, 1e-12);
            TS_ASSERT_DELTA(p_crypt_projection_force->GetB(), 2.0, 1e-12);
            TS_ASSERT_DELTA(p_crypt_projection_force->GetWntChemotaxisStrength(), 15.0, 1e-6);

            delete p_crypt_projection_force;
        }
    }

    void TestForceCollection() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            if (i==60)
            {
                CellPtr p_cell(new Cell(p_apc2, p_model));
                p_cell->SetBirthTime(-10);
                cells.push_back(p_cell);
            }
            else
            {
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetBirthTime(-10);
                cells.push_back(p_cell);
            }
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create two different force laws and add to a std::list
        GeneralisedLinearSpringForce<2> linear_force;

        WntConcentration<2>::Instance()->SetCryptProjectionParameterA(0.0001);
        WntConcentration<2>::Instance()->SetCryptProjectionParameterB(0.0001);
        CryptProjectionForce crypt_projection_force;

        std::vector<AbstractForce<2>* > forces;
        forces.push_back(&linear_force);
        forces.push_back(&crypt_projection_force);

        // Test node force calculation

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }

        // Add force contributions
        for (std::vector<AbstractForce<2>* >::iterator iter = forces.begin();
             iter != forces.end();
             ++iter)
        {
             (*iter)->AddForceContribution(node_forces, cell_population);
        }

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

            TS_ASSERT_DELTA(node_forces[node_index][0], 0.0, 1e-4);
            TS_ASSERT_DELTA(node_forces[node_index][1], 0.0, 1e-4);
        }

        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point = p_mesh->GetNode(59)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];

        p_mesh->SetNode(59, new_point, false);

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > new_node_forces;
        new_node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             new_node_forces.push_back(zero_vector<double>(2));
        }

        // Add force contributions
        for (std::vector<AbstractForce<2>* >::iterator iter = forces.begin();
             iter != forces.end();
             ++iter)
        {
             (*iter)->AddForceContribution(new_node_forces, cell_population);
        }

        // Forces should be twice the forces found using Meineke alone (since a flat crypt is used)
        TS_ASSERT_DELTA(new_node_forces[60][0], 2*0.5*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[60][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[59][0], 2*(-3+4.0/sqrt(7))*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[59][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[58][0], 2*0.5*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[58][1], 0.0, 1e-4);

        WntConcentration<2>::Destroy();
    }

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
        NagaiHondaDifferentialAdhesionForce<2> force; // the true says to use Differential Adhesion.
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

            // Change the value of a CellBasedConfig member
            TS_ASSERT_DELTA(CellBasedConfig::Instance()->GetTransitCellG1Duration(), 2.0, 1e-6);
            CellBasedConfig::Instance()->SetTransitCellG1Duration(15.3);
            TS_ASSERT_DELTA(CellBasedConfig::Instance()->GetTransitCellG1Duration(), 15.3, 1e-6);

            // Serialize via pointer
            NagaiHondaDifferentialAdhesionForce<2>* const p_force = &force;
            output_arch << p_force;

            // Tidy up
            CellBasedConfig::Instance()->Reset();
        }

        {
            // Check CellBasedConfig is reset
            TS_ASSERT_DELTA(CellBasedConfig::Instance()->GetTransitCellG1Duration(), 2.0, 1e-6);

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            NagaiHondaDifferentialAdhesionForce<2>* p_force;

            // Restore from the archive
            input_arch >> p_force;

            // Check CellBasedConfig has been correctly archived
            TS_ASSERT_DELTA(CellBasedConfig::Instance()->GetTransitCellG1Duration(), 15.3, 1e-6);

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

#endif /*TESTFORCESNOTFORRELEASE_HPP_*/

