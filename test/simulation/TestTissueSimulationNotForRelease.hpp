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
#ifndef TESTTISSUESIMULATIONNOTFORRELEASE_HPP_
#define TESTTISSUESIMULATIONNOTFORRELEASE_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>

#include "TissueSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CryptProjectionForce.hpp"
#include "MeshBasedTissueWithGhostNodes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"

/**
 *  Note: Most tests of TissueSimulation are in TestCryptSimulation2d
 */
class TestTissueSimulationNotForRelease : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    /**
     *  Test a tissue simulation with a non-Meineke spring system.
     *
     *  This test consists of a standard crypt projection model simulation with a
     *  radial sloughing cell killer, a crypt projection cell cycle model that
     *  depends on a radial Wnt gradient, and the crypt projection model spring
     *  system, and store the results for use in later archiving tests.
     */
    void TestTissueSimulationWithCryptProjectionSpringSystem() throw (Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();

        double a = 0.2;
        double b = 2.0;
        p_params->SetCryptProjectionParameterA(a);
        p_params->SetCryptProjectionParameterB(b);

        // Set up mesh
        int num_cells_depth = 20;
        int num_cells_width = 20;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        ChasteCuboid<2> bounding_box=p_mesh->CalculateBoundingBox();
        double width_of_mesh = (num_cells_width/(num_cells_width+2.0*thickness_of_ghost_layer))*(bounding_box.GetWidth(0));
        double height_of_mesh = (num_cells_depth/(num_cells_depth+2.0*thickness_of_ghost_layer))*(bounding_box.GetWidth(1));

        p_mesh->Translate(-width_of_mesh/2, -height_of_mesh/2);

        // To start off with, set up all cells to be of type TRANSIT
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            SimpleWntCellCycleModel* p_model = new SimpleWntCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(TRANSIT);
            p_model->SetWntStemThreshold(0.95);

            TissueCellPtr p_cell(new TissueCell(p_state, p_model));
            p_cell->InitialiseCellCycleModel();

            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                  ( p_params->GetTransitCellG1Duration()
                                   +p_params->GetSG2MDuration());
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Make a tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Set up the Wnt gradient
        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Create the force law and pass in to a std::list
        CryptProjectionForce crypt_projection_force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&crypt_projection_force);

        // Make a tissue simulation
        TissueSimulation<2> crypt_projection_simulator(crypt, force_collection, false, false);

        // Create a radial cell killer and pass it in to the tissue simulation
        c_vector<double,2> centre = zero_vector<double>(2);
        double crypt_radius = pow(TissueConfig::Instance()->GetCryptLength()/a, 1.0/b);

        RadialSloughingCellKiller killer(&crypt, centre, crypt_radius);
        crypt_projection_simulator.AddCellKiller(&killer);

        // Set up the simulation
        crypt_projection_simulator.SetOutputDirectory("CryptProjectionSimulation");
        crypt_projection_simulator.SetEndTime(0.25);

        // Run the simulation
        TS_ASSERT_THROWS_NOTHING(crypt_projection_simulator.Solve());

        // These cells just divided and have been gradually moving apart.
        // These results are from time 0.25.
        std::vector<double> node_329_location = crypt_projection_simulator.GetNodeLocation(329);
        std::vector<double> node_494_location = crypt_projection_simulator.GetNodeLocation(494);
        c_vector<double, 2> distance_between;
        distance_between(0) = node_494_location[0] - node_329_location[0];
        distance_between(1) = node_494_location[1] - node_329_location[1];
        TS_ASSERT_DELTA(norm_2(distance_between), 0.6145, 1e-3);

        // Test the Wnt concentration result
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(329)), 0.8753, 1e-3);
        TS_ASSERT_DELTA(p_wnt->GetWntLevel(crypt.GetCellUsingLocationIndex(494)), 0.9175, 1e-3);

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    /**
     * The purpose of this test is to check that it is possible to construct and run
     * a short 1D tissue simulation without throwing any exceptions.
     */
    void Test1dTissueSimulation() throw (Exception)
    {
        // Create mesh
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells so that cell 10 divides at time t=0.5, cell 9 at time t=1.5, etc
        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            TissueCellPtr p_cell(new TissueCell(p_state, p_model));

            double birth_time = -13.5 - i;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Create a tissue
        MeshBasedTissue<1> tissue(mesh, cells);

        // Coverage
        tissue.SetOutputCellIdData(true);

        // Create a force law (no need for a cutoff as we're in 1D)
        GeneralisedLinearSpringForce<1> linear_force;
        std::vector<AbstractForce<1>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up tissue simulation
        TissueSimulation<1> simulator(tissue, force_collection);
        simulator.SetOutputDirectory("Test1DTissueSimulation");
        simulator.SetEndTime(0.6);

        unsigned initial_num_cells = simulator.rGetTissue().GetNumRealCells();
        unsigned initial_num_nodes = simulator.rGetTissue().GetNumNodes();
        unsigned initial_num_elements = (static_cast<MeshBasedTissue<1>* >(&(simulator.rGetTissue())))->rGetMesh().GetNumElements();

        // Run simulation for a short time
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumRealCells(), initial_num_cells + 1);
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumNodes(), initial_num_nodes + 1);
        TS_ASSERT_EQUALS((static_cast<MeshBasedTissue<1>* >(&(simulator.rGetTissue())))->rGetMesh().GetNumElements(), initial_num_elements + 1);
    }

};
#endif /*TESTTISSUESIMULATIONNOTFORRELEASE_HPP_*/
