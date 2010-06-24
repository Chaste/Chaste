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
#include "MeshBasedTissueWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombMutableVertexMeshGenerator.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "ChemotacticForce.hpp"
#include "CellwiseDataGradient.hpp"
#include "CryptProjectionForce.hpp"
#include "NagaiHondaForce.hpp"
#include "WelikyOsterForce.hpp"
#include "VertexCryptBoundaryForce.hpp"
#include "VertexBasedTissue.hpp"
#include "WntConcentration.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "LabelledCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

class TestForcesNotForRelease : public AbstractCellBasedTestSuite
{
public:

    void TestChemotacticForceMethods() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 0;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new LabelledCellMutationState);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCell cell(p_state, p_model);
            cell.SetBirthTime(-10);

            cells.push_back(cell);
        }

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);

        // Set up cellwise data and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(tissue.GetNumRealCells(), 1);
        p_data->SetTissue(&tissue);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            p_data->SetValue(x/50.0, p_mesh->GetNode(i)->GetIndex());
        }

        ChemotacticForce<2> chemotactic_force;

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }
        chemotactic_force.AddForceContribution(node_forces, tissue);

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            unsigned index = tissue.GetLocationIndexUsingCell(*cell_iter);
            double x = tissue.GetLocationOfCellCentre(*cell_iter)[0];
            double c = x/50;
            double norm_grad_c = 1.0/50.0;
            double force_magnitude = chemotactic_force.GetChemotacticForceMagnitude(c, norm_grad_c);

            // Fc = force_magnitude*(1,0), Fspring = 0
            TS_ASSERT_DELTA(node_forces[index][0], force_magnitude, 1e-4);
            TS_ASSERT_DELTA(node_forces[index][1], 0.0, 1e-4);
        }
    }

    void TestChemotacticForceArchiving() throw (Exception)
    {
        // Set up
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "chemotaxis_spring_system.arch";

        {
            // Create mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // SimulationTime is usually set up by a TissueSimulation
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // Create cells
            std::vector<TissueCell> cells;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCell cell(p_state, p_model);
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                cell.SetBirthTime(-50.0);
                cells.push_back(cell);
            }

            // Create tissue
            MeshBasedTissue<2> tissue(mesh,cells);

            // Create force
            ChemotacticForce<2> chemotactic_force;

            // Change the value of a TissueConfig member
            TS_ASSERT_DELTA(TissueConfig::Instance()->GetStemCellG1Duration(), 14.0, 1e-6);
            TissueConfig::Instance()->SetStemCellG1Duration(17.68);
            TS_ASSERT_DELTA(TissueConfig::Instance()->GetStemCellG1Duration(), 17.68, 1e-6);

            // Serialize force (and hence TissueConfig) via pointer
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            ChemotacticForce<2>* const p_chemotactic_force = &chemotactic_force;
            output_arch << p_chemotactic_force;

            // Tidy up
            TissueConfig::Instance()->Reset();
        }

        {
            // Check TissueConfig is reset
            TS_ASSERT_DELTA(TissueConfig::Instance()->GetStemCellG1Duration(), 14.0, 1e-6);

            ArchiveLocationInfo::SetMeshPathname("mesh/test/data/", "square_2_elements");

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore force (and hence TissueConfig) from the archive
            ChemotacticForce<2>* p_chemotactic_force;
            input_arch >> p_chemotactic_force;

            // Check TissueConfig has been correctly archived
            TS_ASSERT_DELTA(TissueConfig::Instance()->GetStemCellG1Duration(), 17.68, 1e-6);

            // Tidy up
            delete p_chemotactic_force;
        }
    }

    void TestCryptProjectionForceMethods() throw (Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

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
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCell cell(p_state, p_model);

            if (i==4 || i==5)
            {
                cell.SetBirthTime(-0.5);
            }
            else
            {
                cell.SetBirthTime(-10.0);
            }
            cells.push_back(cell);
        }

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        std::set<TissueCell*> cell_pair_4_5 = tissue.CreateCellPair(tissue.rGetCellUsingLocationIndex(4), tissue.rGetCellUsingLocationIndex(5));
        tissue.MarkSpring(cell_pair_4_5);

        // Create a spring system with crypt surface z = 2*r
        p_params->SetCryptProjectionParameterA(2.0);
        p_params->SetCryptProjectionParameterB(1.0);
        CryptProjectionForce crypt_projection_force;

        // Test get methods
        TS_ASSERT_DELTA(crypt_projection_force.GetA(), 2.0, 1e-12);
        TS_ASSERT_DELTA(crypt_projection_force.GetB(), 1.0, 1e-12);

        // Test crypt height and gradient calculations
        c_vector<double, 2> node_location_2d = p_mesh->GetNode(0)->rGetLocation();
        TS_ASSERT_DELTA(crypt_projection_force.CalculateCryptSurfaceHeightAtPoint(node_location_2d), 2.0*pow(norm_2(node_location_2d),1.0), 1e-12);
        TS_ASSERT_DELTA(crypt_projection_force.CalculateCryptSurfaceDerivativeAtPoint(node_location_2d), 2.0, 1e-12);

        // Test updating of mNode3dLocationMap
        crypt_projection_force.UpdateNode3dLocationMap(tissue);

        // Move a node slightly
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = node_location_2d[0]+0.05;
        new_point.rGetLocation()[1] = node_location_2d[1];
        p_mesh->SetNode(0, new_point, false);


        // Test UpdateNode3dLocationMap()

        c_vector<double, 2> new_node_location_2d;
        new_node_location_2d[0] = new_point.rGetLocation()[0];
        new_node_location_2d[1] = new_point.rGetLocation()[1];

        crypt_projection_force.UpdateNode3dLocationMap(tissue);

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
                                                                            tissue);

        TS_ASSERT_DELTA(force_on_spring[0], -5.7594, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1],  0.0230, 1e-4);


        // Test force calculation with a cutoff

        double dist = norm_2( p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(),
                              p_element->GetNode(1)->rGetLocation()) );

        crypt_projection_force.UseCutoffPoint(dist - 0.1);

        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                            p_element->GetNodeGlobalIndex(0),
                                                                            tissue);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);


        // Test force calculation for a pair of newly born neighbouring cells
        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(4, 5, tissue);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);

        tissue.UnmarkSpring(cell_pair_4_5);

        // For coverage, test force calculation for a pair of neighbouring apoptotic cells
        tissue.rGetCellUsingLocationIndex(6).StartApoptosis();
        tissue.rGetCellUsingLocationIndex(7).StartApoptosis();
        force_on_spring = crypt_projection_force.CalculateForceBetweenNodes(6, 7, tissue);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);

        // Test force calculation for a particular node

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }

        crypt_projection_force.AddForceContribution(node_forces, tissue);

        TS_ASSERT_DELTA(node_forces[0][0], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_forces[0][1], 0.0, 1e-4);

        // Test that in the case of a flat crypt surface (mA=mB=0), the results are the same as for Meineke2001SpringSystem
        p_params->SetCryptProjectionParameterA(0.001);
        p_params->SetCryptProjectionParameterB(0.001);
        CryptProjectionForce flat_crypt_projection_force;
        GeneralisedLinearSpringForce<2> linear_force;

        // Normally this would be set up at the start of rCalculateforcesOfEachNode
        flat_crypt_projection_force.UpdateNode3dLocationMap(tissue);

        for (MeshBasedTissue<2>::SpringIterator spring_iterator = tissue.SpringsBegin();
            spring_iterator != tissue.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

            c_vector<double, 2> force_flat = flat_crypt_projection_force.CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, tissue);
            c_vector<double, 2> force_meineke = linear_force.CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, tissue);

            TS_ASSERT_DELTA( force_flat[0], force_meineke[0], 1e-3);
            TS_ASSERT_DELTA( force_flat[1], force_meineke[1], 1e-3);
        }
    }

    /**
     * Note: WntBasedChemotaxis should be possible in other force laws. If/when
     * this is implemented, this test should be moved to somewhere more appropriate.
     */
    void TestCryptProjectionForceWithWntBasedChemotaxis() throw (Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

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
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCell cell(p_state, p_model);
            cell.SetBirthTime(-10.0);

            cells.push_back(cell);
        }

        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        std::set<TissueCell*> cell_pair_4_5 = tissue.CreateCellPair(tissue.rGetCellUsingLocationIndex(4), tissue.rGetCellUsingLocationIndex(5));
        tissue.MarkSpring(cell_pair_4_5);

        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetTissue(tissue);

        // Create a spring system with crypt surface z = 2*r
        p_params->SetCryptProjectionParameterA(2.0);
        p_params->SetCryptProjectionParameterB(1.0);
        CryptProjectionForce crypt_projection_force;

        crypt_projection_force.SetWntChemotaxis(false);

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > old_node_forces;
        old_node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             old_node_forces.push_back(zero_vector<double>(2));
        }

        // Calculate node forces
        crypt_projection_force.AddForceContribution(old_node_forces, tissue);

        // Store the force of a particular node without Wnt-chemotaxis
        c_vector<double,2> old_force = old_node_forces[11];

        // Now turn on Wnt-chemotaxis
        crypt_projection_force.SetWntChemotaxis(true);

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > new_node_forces;
        new_node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             new_node_forces.push_back(zero_vector<double>(2));
        }

        // Calculate node forces
        crypt_projection_force.AddForceContribution(new_node_forces, tissue);

        // Store the force of the same node, but now with Wnt-chemotaxis
        c_vector<double,2> new_force = new_node_forces[11];

        double wnt_chemotaxis_strength = TissueConfig::Instance()->GetWntChemotaxisStrength();
        c_vector<double,2> wnt_component = wnt_chemotaxis_strength*WntConcentration<2>::Instance()->GetWntGradient(cells[11]);

        TS_ASSERT_DELTA(new_force[0], old_force[0]+wnt_component[0], 1e-4);
        TS_ASSERT_DELTA(new_force[1], old_force[1]+wnt_component[1], 1e-4);
    }

    void TestCryptProjectionForceWithArchiving() throw (Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "crypt_projection_spring_system.arch";

        {
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");

            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

            std::vector<TissueCell> cells;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCell cell(p_state, p_model);
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                cell.SetBirthTime(-50.0);
                cells.push_back(cell);
            }

            MeshBasedTissue<2> crypt(mesh,cells);
            p_params->SetCryptProjectionParameterA(1.0);
            p_params->SetCryptProjectionParameterB(2.0);
            CryptProjectionForce crypt_projection_force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            CryptProjectionForce* const p_crypt_projection_force = &crypt_projection_force;

            p_crypt_projection_force->UseCutoffPoint(1.1);

            output_arch << p_crypt_projection_force;
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

            delete p_crypt_projection_force;
        }
    }

    void TestForceCollection() throw (Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            if (i==60)
            {
                TissueCell cell(p_apc2, p_model);
                cell.SetBirthTime(-10);
                cells.push_back(cell);
            }
            else
            {
                TissueCell cell(p_state, p_model);
                cell.SetBirthTime(-10);
                cells.push_back(cell);
            }
        }

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);

        // Create two different force laws and add to a std::list
        GeneralisedLinearSpringForce<2> linear_force;

        p_params->SetCryptProjectionParameterA(0.0001);
        p_params->SetCryptProjectionParameterB(0.0001);
        CryptProjectionForce crypt_projection_force;

        std::vector<AbstractForce<2>* > forces;
        forces.push_back(&linear_force);
        forces.push_back(&crypt_projection_force);

        // Test node force calculation

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }

        // Add force contributions
        for (std::vector<AbstractForce<2>* >::iterator iter = forces.begin();
             iter != forces.end();
             ++iter)
        {
             (*iter)->AddForceContribution(node_forces, tissue);
        }

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            unsigned node_index = tissue.GetLocationIndexUsingCell(*cell_iter);

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
        new_node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
             new_node_forces.push_back(zero_vector<double>(2));
        }

        // Add force contributions
        for (std::vector<AbstractForce<2>* >::iterator iter = forces.begin();
             iter != forces.end();
             ++iter)
        {
             (*iter)->AddForceContribution(new_node_forces, tissue);
        }

        // Forces should be twice the forces found using Meineke alone (since a flat crypt is used)
        TS_ASSERT_DELTA(new_node_forces[60][0], 2*0.5*p_params->GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[60][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[59][0], 2*(-3+4.0/sqrt(7))*p_params->GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[59][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[58][0], 2*0.5*p_params->GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[58][1], 0.0, 1e-4);
    }

    void TestNagaiHondaForceMethods() throw (Exception)
    {
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        std::vector<double> angles = std::vector<double>(num_nodes);

        for (unsigned i=0; i<num_nodes; i++)
        {
            angles[i] = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(angles[i]), sin(angles[i])));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Set up the cell
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(DIFFERENTIATED);

        TissueCell cell(p_state, p_model);
        cell.SetBirthTime(-1.0);
        cells.push_back(cell);

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);
        tissue.InitialiseCells();

        // Create a force system
        NagaiHondaForce<2> force;

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        force.AddForceContribution(node_forces, tissue);

        // The force on each node should be radially inward, with the same magnitude for all nodes
        double force_magnitude = norm_2(node_forces[0]);

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(node_forces[i]), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(node_forces[i][0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(node_forces[i][1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        double normal_target_area = tissue.GetTargetAreaOfCell(tissue.rGetCellUsingLocationIndex(0));

        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.25, 2);

        // Set the cell to be necrotic
        tissue.rGetCellUsingLocationIndex(0).StartApoptosis();

        double initial_apoptotic_target_area = tissue.GetTargetAreaOfCell(tissue.rGetCellUsingLocationIndex(0));

        TS_ASSERT_DELTA(normal_target_area, initial_apoptotic_target_area, 1e-6);

        // Reset force vector
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            node_forces[i] = zero_vector<double>(2);
        }

        force.AddForceContribution(node_forces, tissue);

        // The force on each node should not yet be affected by setting the cell to be apoptotic
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(node_forces[i]), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(node_forces[i][0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(node_forces[i][1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Increment time
        p_simulation_time->IncrementTimeOneStep();

        double later_apoptotic_target_area = tissue.GetTargetAreaOfCell(tissue.rGetCellUsingLocationIndex(0));

        TS_ASSERT_LESS_THAN(later_apoptotic_target_area, initial_apoptotic_target_area);

        TS_ASSERT_DELTA(tissue.rGetCellUsingLocationIndex(0).GetTimeUntilDeath(), 0.125, 1e-6);

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            node_forces[i] = zero_vector<double>(2);
        }

        force.AddForceContribution(node_forces, tissue);

        // Now the forces should be affected
        double apoptotic_force_magnitude = norm_2(node_forces[0]);
        TS_ASSERT_LESS_THAN(force_magnitude, apoptotic_force_magnitude);
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(node_forces[i]), apoptotic_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(node_forces[i][0], -apoptotic_force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(node_forces[i][1], -apoptotic_force_magnitude*sin(angles[i]), 1e-4);
        }
    }

    void TestWelikyOsterForceMethods() throw (Exception)
    {
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        std::vector<double> angles = std::vector<double>(num_nodes);

        for (unsigned i=0; i<num_nodes; i++)
        {
            angles[i] = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(angles[i]), sin(angles[i])));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Set up the cell
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(DIFFERENTIATED);

        TissueCell cell(p_state, p_model);
        cell.SetBirthTime(-1.0);
        cells.push_back(cell);

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);
        tissue.InitialiseCells();

        // Create a force system
        WelikyOsterForce<2> force;

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        force.AddForceContribution(node_forces, tissue);

        // The force on each node should be radially inward, with the same magnitude for all nodes
        double force_magnitude = norm_2(node_forces[0]);

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(node_forces[i]), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(node_forces[i][0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(node_forces[i][1], -force_magnitude*sin(angles[i]), 1e-4);
        }
    }

    void TestArchivingWelikyOsterForce() throw (Exception)
    {
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "weliky_oster.arch";

        {
            // Construct a 2D vertex mesh consisting of a single element
            std::vector<Node<2>*> nodes;
            unsigned num_nodes = 20;
            std::vector<double> angles = std::vector<double>(num_nodes);

            for (unsigned i=0; i<num_nodes; i++)
            {
                angles[i] = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
                nodes.push_back(new Node<2>(i, false, cos(angles[i]), sin(angles[i])));
            }

            std::vector<VertexElement<2,2>*> elements;
            elements.push_back(new VertexElement<2,2>(0, nodes));

            double cell_swap_threshold = 0.01;
            double edge_division_threshold = 2.0;
            MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

            // Set up the cell
            std::vector<TissueCell> cells;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            TissueCell cell(p_state, p_model);
            cell.SetBirthTime(-1.0);
            cells.push_back(cell);

            // Create tissue
            VertexBasedTissue<2> tissue(mesh, cells);
            tissue.InitialiseCells();

            // Create a force system
            WelikyOsterForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Change the value of a TissueConfig member
            TS_ASSERT_DELTA(TissueConfig::Instance()->GetTransitCellG1Duration(), 2.0, 1e-6);
            TissueConfig::Instance()->SetTransitCellG1Duration(15.3);
            TS_ASSERT_DELTA(TissueConfig::Instance()->GetTransitCellG1Duration(), 15.3, 1e-6);

            // Serialize via pointer
            WelikyOsterForce<2>* const p_force = &force;
            output_arch << p_force;

            // Tidy up
            TissueConfig::Instance()->Reset();
        }

        {
            // Check TissueConfig is reset
            TS_ASSERT_DELTA(TissueConfig::Instance()->GetTransitCellG1Duration(), 2.0, 1e-6);

//            ArchiveLocationInfo::SetMeshPathname("mesh/test/data/", "square_2_elements");

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            WelikyOsterForce<2>* p_force;

            // Restore from the archive
            input_arch >> p_force;

            // Check TissueConfig has been correctly archived
            TS_ASSERT_DELTA(TissueConfig::Instance()->GetTransitCellG1Duration(), 15.3, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestVertexCryptBoundaryForce() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        HoneycombMutableVertexMeshGenerator generator(5, 5, false, 0.1, 0.5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMutableMesh();

        // Translate mesh so that some points are below y=0
        p_mesh->Translate(0.0, -3.0);

        // Set up cells, one for each VertexElement. Give each cell
        // a birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            TissueCell cell(p_state, p_model);
            double birth_time = 0.0 - elem_index;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(*p_mesh, cells);

        // Create a force system
        VertexCryptBoundaryForce<2> force(100);

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        force.AddForceContribution(node_forces, tissue);

        // Check forces are correct
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(node_forces[i][0], 0.0, 1e-4);

            double y = tissue.GetNode(i)->rGetLocation()[1];
            if (y >= 0.0)
            {
                // If y > 0, the force contribution should be zero...
                TS_ASSERT_DELTA(node_forces[i][1], 0.0, 1e-4);
            }
            else
            {
                // ...otherwise, the force contribution should be quadratic in y
                double expected_force = force.GetForceStrength()*y*y;
                TS_ASSERT_DELTA(node_forces[i][1], expected_force, 1e-4);
            }
        }
    }

};

#endif /*TESTFORCESNOTFORRELEASE_HPP_*/

