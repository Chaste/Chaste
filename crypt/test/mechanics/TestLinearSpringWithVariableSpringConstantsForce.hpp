/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTLINEARSPRINGWITHVARIABLESPRINGCONSTANTSFORCE_HPP_
#define TESTLINEARSPRINGWITHVARIABLESPRINGCONSTANTSFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CryptCellsGenerator.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
//#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "WntConcentration.hpp"
#include "CryptSimulation2d.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

class TestLinearSpringWithVariableSpringConstantsForce : public AbstractCellBasedTestSuite
{
public:

    void TestGeneralisedLinearSpringForceWithSpringConstantsForIngeBCatCellsSimulationNoGhostsFails()
    {
        double crypt_length = 1.1*12.0*sqrt(3)/2.0;
      //This one is going to throw because the ghost layers haven't been set up properly
        CylindricalHoneycombMeshGenerator generator(6, 12, 0, 1.1);

        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, false);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Set up simulation
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("TestLinearSpringWithVariableSpringConstantsForceSimulation");
        simulator.SetEndTime(0.1);

        MAKE_PTR(LinearSpringWithVariableSpringConstantsForce<2>, p_linear_force);
        p_linear_force->SetBetaCateninSprings(true);
        crypt.CreateVoronoiTessellation();  // this method is normally called in a simulation loop

        simulator.AddForce(p_linear_force);

        TS_ASSERT_THROWS_THIS(simulator.Solve(), "Spring iterator tried to calculate interaction between degenerate cells on the boundary of the mesh.  Have you set ghost layers correctly?");

        // Tidy up
        WntConcentration<2>::Destroy();
        SimulationTime::Destroy();
    }
    void TestGeneralisedLinearSpringForceWithSpringConstantsForIngeBCatCellsSimulationWithGhosts()
    {
        double crypt_length = 1.1*12.0*sqrt(3)/2.0;
        //Do it again, but with ghost layers set up properly
        CylindricalHoneycombMeshGenerator generator(6, 12, 2, 1.1);
        //                                                 ^-- That was a zero before
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, false);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Set up simulation
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("TestLinearSpringWithVariableSpringConstantsForceSimulation");
        simulator.SetEndTime(0.1);

        MAKE_PTR(LinearSpringWithVariableSpringConstantsForce<2>, p_linear_force);

        p_linear_force->SetBetaCateninSprings(true);
        crypt.CreateVoronoiTessellation();  // this method is normally called in a simulation loop

        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Tidy up
        WntConcentration<2>::Destroy();
        SimulationTime::Destroy();
    }


    void TestGeneralisedLinearSpringForceWithSpringConstantsForIngeBCatCells()
    {
        double crypt_length = 1.1*12.0*sqrt(3)/2.0;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        CylindricalHoneycombMeshGenerator generator(6, 12, 0, 1.1);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, false);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // As there is no cell-based simulation, we must explicitly initialise the cells
        crypt.InitialiseCells();

        LinearSpringWithVariableSpringConstantsForce<2> linear_force;

        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.50, 1e-10);

        linear_force.SetBetaCateninSprings(true);
        crypt.CreateVoronoiTessellation();  // this method is normally called in a simulation loop

        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.5*8.59312/18.14, 1e-5);

        linear_force.SetBetaCatSpringScaler(20/6.0);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(20, 21, crypt)), 1.5*8.59312/20.0, 1e-5);

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    void TestGeneralisedLinearSpringForceWithEdgeLengthBasedSpring() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, crypt_width/cells_across);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CellPropertyRegistry::Instance()->Clear();
        RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

        // Set up cells
        std::vector<CellPtr> cells;
        cells.clear();
        unsigned num_cells = location_indices.empty() ? p_mesh->GetNumNodes() : location_indices.size();
        cells.reserve(num_cells);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellProliferativeType cell_type;
            unsigned generation;
            double y = 0.0;

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            }

            FixedDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel;
            p_cell_cycle_model->SetDimension(2);

            double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
            double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

            double birth_time = 0.0;
            birth_time = -p_random_num_gen->ranf();

            if (y <= 0.3)
            {
                cell_type = STEM;
                generation = 0;
                birth_time *= typical_stem_cycle_time; // hours
            }
            else if (y < 2.0)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 3.0)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 4.0)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else
            {
                cell_type = p_cell_cycle_model->CanCellTerminallyDifferentiate() ? DIFFERENTIATED : TRANSIT;
                generation = 4;
                birth_time *= typical_transit_cycle_time; // hours
            }

            p_cell_cycle_model->SetGeneration(generation);
            p_cell_cycle_model->SetCellProliferativeType(cell_type);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetBirthTime(birth_time);

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                cells.push_back(p_cell);
            }
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        LinearSpringWithVariableSpringConstantsForce<2> linear_force;

        // Check that the force between nodes is correctly calculated when the 'spring constant' is constant
        linear_force.SetEdgeBasedSpringConstant(false);

        for (MeshBasedCellPopulation<2>::SpringIterator spring_iterator = cell_population.SpringsBegin();
            spring_iterator != cell_population.SpringsEnd();
            ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                nodeB_global_index,
                                                                                cell_population);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 6.25, 1e-3);
        }

        // Check that the force between nodes is correctly calculated when the 'spring constant'
        // is proportional to the length of the edge between adjacent cells
        linear_force.SetEdgeBasedSpringConstant(true);
        cell_population.CreateVoronoiTessellation();  // normally done in a simulation loop

        for (MeshBasedCellPopulation<2>::SpringIterator spring_iterator = cell_population.SpringsBegin();
             spring_iterator != cell_population.SpringsEnd();
             ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                nodeB_global_index,
                                                                                cell_population);

            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 4.34027778, 1e-3);
        }

        // Choose two interior neighbour nodes
        c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(41, 42, cell_population);
        TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 4.34027778, 1e-3);

        // Now move node 42 a bit and check that the force calculation changes correctly
        c_vector<double,2> shift;
        shift[0] = 0.1;
        shift[1] = 0.0;
        ChastePoint<2> new_point(p_mesh->GetNode(42u)->rGetLocation() + shift);
        p_mesh->SetNode(21, new_point, false);

        // Check that the new force between nodes is correctly calculated
        cell_population.CreateVoronoiTessellation();
        c_vector<double, 2> new_force = linear_force.CalculateForceBetweenNodes(41, 42, cell_population);

        // Force calculation: shift is along x-axis so we should have
        // new_edge_length = (5/6 + shift[0])*tan(0.5*arctan(5*sqrt(3)/(5 + 12*shift[0]))),
        // force^2 = mu^2 * (new_edge_length*sqrt(3))^2 * (1 - 5/6 - shift[0])^2
        TS_ASSERT_DELTA(new_force[0]*new_force[0] + new_force[1]*new_force[1], 4.34024, 1e-3);
    }

    void TestGeneralisedLinearSpringForceWithEdgeBasedSpringsOnPeriodicMesh() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 3;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CellPropertyRegistry::Instance()->Clear();
        RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

        // Set up cells
        std::vector<CellPtr> cells;
        cells.clear();
        unsigned num_cells = location_indices.empty() ? p_mesh->GetNumNodes() : location_indices.size();
        cells.reserve(num_cells);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellProliferativeType cell_type;
            unsigned generation;
            double y = 0.0;

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            }

            FixedDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel;
            p_cell_cycle_model->SetDimension(2);

            double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
            double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

            double birth_time = 0.0;
            birth_time = -p_random_num_gen->ranf();

            if (y <= 0.3)
            {
                cell_type = STEM;
                generation = 0;
                birth_time *= typical_stem_cycle_time; // hours
            }
            else if (y < 2.0)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 3.0)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 4.0)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else
            {
                cell_type = p_cell_cycle_model->CanCellTerminallyDifferentiate() ? DIFFERENTIATED : TRANSIT;
                generation = 4;
                birth_time *= typical_transit_cycle_time; // hours
            }

            p_cell_cycle_model->SetGeneration(generation);
            p_cell_cycle_model->SetCellProliferativeType(cell_type);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetBirthTime(birth_time);

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                cells.push_back(p_cell);
            }
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        LinearSpringWithVariableSpringConstantsForce<2> linear_force;

        // Check that the force between nodes is correctly calculated when the spring constant is constant (!)
        linear_force.SetEdgeBasedSpringConstant(false);

        for (MeshBasedCellPopulation<2>::SpringIterator spring_iterator = cell_population.SpringsBegin();
             spring_iterator != cell_population.SpringsEnd();
             ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                nodeB_global_index,
                                                                                cell_population);

            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 6.25, 1e-3);
        }

        // Check that the force between nodes is correctly calculated when the spring constant
        // is proportional to the length of the edge between adjacenet cells
        linear_force.SetEdgeBasedSpringConstant(true);
        cell_population.CreateVoronoiTessellation();

        for (MeshBasedCellPopulation<2>::SpringIterator spring_iterator = cell_population.SpringsBegin();
             spring_iterator != cell_population.SpringsEnd();
             ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
            c_vector<double, 2> force = linear_force.CalculateForceBetweenNodes(nodeA_global_index,
                                                                                nodeB_global_index,
                                                                                cell_population);
            TS_ASSERT_DELTA(force[0]*force[0] + force[1]*force[1], 4.34027778, 1e-3);
        }
    }

    void TestGeneralisedLinearSpringForceWithSpringConstantsForMutantCells()
    {
        // Create a small cell population
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0));
        nodes.push_back(new Node<2>(1, false, 0, 2));
        nodes.push_back(new Node<2>(2, false, 2, 2));
        nodes.push_back(new Node<2>(3, false, 2, 0));

        MutableMesh<2,2> mesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        LinearSpringWithVariableSpringConstantsForce<2> linear_force;

        // Set cells' mutation states
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_apc2(new ApcTwoHitCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_bcat1(new BetaCateninOneHitCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);

        cell_population.GetCellUsingLocationIndex(0)->SetMutationState(p_state);
        cell_population.GetCellUsingLocationIndex(1)->AddCellProperty(p_label);
        cell_population.GetCellUsingLocationIndex(2)->SetMutationState(p_apc2);
        cell_population.GetCellUsingLocationIndex(3)->SetMutationState(p_bcat1);

        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(0, 1, cell_population)), 15.0, 1e-10);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(1, 2, cell_population)), 15.0, 1e-10);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(2, 3, cell_population)), 15.0, 1e-10);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(3, 0, cell_population)), 15.0, 1e-10);

        linear_force.SetMutantSprings(true);

        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(0, 1, cell_population)), 15.0, 1e-10);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(1, 2, cell_population)), 22.5, 1e-10);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(2, 3, cell_population)), 30.0, 1e-10);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(3, 0, cell_population)), 22.5, 1e-10);

        linear_force.SetMutantSprings(true, 4.0, 3.0);

        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(0, 1, cell_population)), 15.0, 1e-10);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(1, 2, cell_population)), 45.0, 1e-10);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(2, 3, cell_population)), 60.0, 1e-10);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(3, 0, cell_population)), 45.0, 1e-10);
    }

    void TestGeneralisedLinearSpringForceWithSpringConstantsForApoptoticCells()
    {
        // Set up stretched cell population
        HoneycombMeshGenerator generator(4, 4, 0, 2.0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        MeshBasedCellPopulationWithGhostNodes<2> stretched_cell_population(*p_mesh, cells, location_indices);

        // As there is no cell-based simulation we must explicitly initialise the cells
        stretched_cell_population.InitialiseCells();

        // Set one of the non-boundary cells to be necrotic
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(new ApoptoticCellProperty);
        stretched_cell_population.GetCellUsingLocationIndex(6)->AddCellProperty(p_apoptotic_state);

        LinearSpringWithVariableSpringConstantsForce<2> linear_force;
        linear_force.SetApoptoticSprings(true);

        TS_ASSERT_EQUALS(stretched_cell_population.GetCellUsingLocationIndex(6)->HasCellProperty<ApoptoticCellProperty>(), true);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(6, 10, stretched_cell_population)), 3.3333, 1e-4);

        // Set a neighbouring cell to be necrotic
        stretched_cell_population.GetCellUsingLocationIndex(10)->AddCellProperty(p_apoptotic_state);

        TS_ASSERT_EQUALS(stretched_cell_population.GetCellUsingLocationIndex(10)->HasCellProperty<ApoptoticCellProperty>(), true);
        TS_ASSERT_DELTA(norm_2(linear_force.CalculateForceBetweenNodes(6, 10, stretched_cell_population)), 1.8750, 1e-4);

        // Now do similar tests for a squashed cell population
        HoneycombMeshGenerator generator2(4, 4, 0, 0.5);
        MutableMesh<2,2>* p_mesh2 = generator2.GetMesh();
        std::vector<unsigned> location_indices2 = generator2.GetCellLocationIndices();

        std::vector<CellPtr> cells2;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel,2> cells_generator2;
        cells_generator2.GenerateGivenLocationIndices(cells2, location_indices2);

        MeshBasedCellPopulationWithGhostNodes<2> squashed_cell_population(*p_mesh2, cells2, location_indices2);
        squashed_cell_population.InitialiseCells();

        squashed_cell_population.GetCellUsingLocationIndex(6)->AddCellProperty(p_apoptotic_state);

        LinearSpringWithVariableSpringConstantsForce<2> linear_force2;

        // Test set/get methods
        TS_ASSERT_DELTA(linear_force2.GetApoptoticSpringTensionStiffness(), 3.75, 1e-6);
        TS_ASSERT_DELTA(linear_force2.GetApoptoticSpringCompressionStiffness(), 11.25, 1e-6);

        linear_force2.SetApoptoticSprings(true);

        TS_ASSERT_DELTA(norm_2(linear_force2.CalculateForceBetweenNodes(6, 10, squashed_cell_population)), 4.0909, 1e-4);

        squashed_cell_population.GetCellUsingLocationIndex(10)->AddCellProperty(p_apoptotic_state);

        TS_ASSERT_DELTA(norm_2(linear_force2.CalculateForceBetweenNodes(6, 10, squashed_cell_population)), 2.8125, 1e-4);
    }

    void TestForceOutputParameters()
    {
        std::string output_directory = "TestForcesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with LinearSpringWithVariableSpringConstantsForce
        LinearSpringWithVariableSpringConstantsForce<2> variable_force;
        variable_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(variable_force.GetIdentifier(), "LinearSpringWithVariableSpringConstantsForce-2");

        out_stream variable_force_parameter_file = output_file_handler.OpenOutputFile("variable_results.parameters");
        variable_force.OutputForceParameters(variable_force_parameter_file);
        variable_force_parameter_file->close();

        std::string variable_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + variable_force_results_dir + "variable_results.parameters  crypt/test/data/TestForcesForCrypt/variable_results.parameters").c_str()), 0);
    }

    void TestLinearSpringWithVariableSpringConstantsForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "meineke_spring_system.arch";

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

            MeshBasedCellPopulation<2> cell_population(mesh, cells);
            LinearSpringWithVariableSpringConstantsForce<2> linear_force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            LinearSpringWithVariableSpringConstantsForce<2>* const p_linear_force = &linear_force;

            p_linear_force->SetCutOffLength(1.1);
            p_linear_force->SetEdgeBasedSpringConstant(true);
            p_linear_force->SetMutantSprings(true, 0.2, 0.3);
            p_linear_force->SetBetaCateninSprings(true);
            p_linear_force->SetApoptoticSprings(true);
            p_linear_force->SetBetaCatSpringScaler(20/6.0);
            p_linear_force->SetApoptoticSpringTensionStiffness(4.0);
            p_linear_force->SetApoptoticSpringCompressionStiffness(5.0);

            output_arch << p_linear_force;
        }

        {
            ArchiveLocationInfo::SetMeshPathname("mesh/test/data", "square_2_elements");

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            LinearSpringWithVariableSpringConstantsForce<2>* p_linear_force;

            // Restore from the archive
            input_arch >> p_linear_force;

            // Test the member data
            TS_ASSERT_EQUALS(p_linear_force->mUseCutOffLength, true);
            TS_ASSERT_EQUALS(p_linear_force->mMechanicsCutOffLength, 1.1);
            TS_ASSERT_EQUALS(p_linear_force->mUseEdgeBasedSpringConstant, true);
            TS_ASSERT_EQUALS(p_linear_force->mUseEdgeBasedSpringConstant, true);
            TS_ASSERT_EQUALS(p_linear_force->mUseMutantSprings, true);
            TS_ASSERT_EQUALS(p_linear_force->mUseBCatSprings, true);
            TS_ASSERT_EQUALS(p_linear_force->mUseApoptoticSprings, true);
            TS_ASSERT_DELTA(p_linear_force->mMutantMutantMultiplier, 0.2, 1e-12);
            TS_ASSERT_DELTA(p_linear_force->mNormalMutantMultiplier, 0.3, 1e-12);
            TS_ASSERT_DELTA(p_linear_force->mBetaCatSpringScaler, 20/6.0, 1e-12);
            TS_ASSERT_DELTA(p_linear_force->mApoptoticSpringTensionStiffness, 4.0, 1e-12);
            TS_ASSERT_DELTA(p_linear_force->mApoptoticSpringCompressionStiffness, 5.0, 1e-12);

            delete p_linear_force;
        }
    }

    void TestLinearSpringWithVariableSpringConstantsForceWithNodeBasedCellPopulation() throw (Exception)
    {
        // Create a NodeBasedCellPopulation
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 10;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double x = (double)(i);
            double y = (double)(i);
            nodes.push_back(new Node<2>(i, true, x, y));
        }

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, num_nodes);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a vector forces on nodes
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(num_nodes);
        for (unsigned i=0; i<num_nodes; i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        // Test that LinearSpringWithVariableSpringConstantsForce throws the correct exception
        LinearSpringWithVariableSpringConstantsForce<2> spring_force;
        TS_ASSERT_THROWS_THIS(spring_force.AddForceContribution(node_forces, cell_population),
                "LinearSpringWithVariableSpringConstantsForce is to be used with a subclass of MeshBasedCellPopulation only");

        // When the node-only mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /*TESTLINEARSPRINGWITHVARIABLESPRINGCONSTANTSFORCE_HPP_*/
