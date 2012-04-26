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

#ifndef TESTFORCES_HPP_
#define TESTFORCES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "GeneralisedLinearSpringForce.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "ChemotacticForce.hpp"
#include "RepulsionForce.hpp"
#include "NagaiHondaForce.hpp"
#include "WelikyOsterForce.hpp"
#include "DiffusionForce.hpp"
#include "CellwiseData.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"

class TestForces : public AbstractCellBasedTestSuite
{
public:

    void TestGeneralisedLinearSpringForceMethods() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size(), location_indices);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create force
        GeneralisedLinearSpringForce<2> linear_force;

        // Test set/get method
        TS_ASSERT_DELTA(linear_force.GetMeinekeDivisionRestingSpringLength(), 0.5, 1e-6);
        TS_ASSERT_DELTA(linear_force.GetMeinekeSpringStiffness(), 15.0, 1e-6);
        TS_ASSERT_DELTA(linear_force.GetMeinekeSpringGrowthDuration(), 1.0, 1e-6);
        TS_ASSERT_EQUALS(linear_force.GetUseCutOffLength(), false);
        TS_ASSERT_DELTA(linear_force.GetCutOffLength(), DBL_MAX, 1e-6);

        linear_force.SetMeinekeDivisionRestingSpringLength(0.8);
        linear_force.SetMeinekeSpringStiffness(20.0);
        linear_force.SetMeinekeSpringGrowthDuration(2.0);
        linear_force.SetCutOffLength(1.5);

        TS_ASSERT_DELTA(linear_force.GetMeinekeDivisionRestingSpringLength(), 0.8, 1e-6);
        TS_ASSERT_DELTA(linear_force.GetMeinekeSpringStiffness(), 20.0, 1e-6);
        TS_ASSERT_DELTA(linear_force.GetMeinekeSpringGrowthDuration(), 2.0, 1e-6);
        TS_ASSERT_EQUALS(linear_force.GetUseCutOffLength(), true);
        TS_ASSERT_DELTA(linear_force.GetCutOffLength(), 1.5, 1e-6);

        linear_force.SetMeinekeDivisionRestingSpringLength(0.5);
        linear_force.SetMeinekeSpringStiffness(15.0);
        linear_force.SetMeinekeSpringGrowthDuration(1.0);

        // Reset cut off length
        linear_force.SetCutOffLength(DBL_MAX);

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }

        // Test node force calculation
        linear_force.AddForceContribution(node_forces, cell_population);

        // Test forces on non-ghost nodes
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
        linear_force.AddForceContribution(new_node_forces, cell_population);

        TS_ASSERT_DELTA(new_node_forces[60][0], 0.5*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[60][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[59][0], (-3+4.0/sqrt(7))*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[59][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(new_node_forces[58][0], 0.5*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(new_node_forces[58][1], 0.0, 1e-4);

        // Test spring force calculation
        c_vector<double,2> force_on_spring; // between nodes 59 and 60

        // Find one of the elements that nodes 59 and 60 live on
        ChastePoint<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0] + 0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01;

        unsigned elem_index = p_mesh->GetContainingElementIndex(new_point2, false);
        Element<2,2>* p_element = p_mesh->GetElement(elem_index);

        force_on_spring = linear_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                  p_element->GetNodeGlobalIndex(0),
                                                                  cell_population);

        TS_ASSERT_DELTA(force_on_spring[0], 0.5*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);

        // Test force with cutoff point
        double dist = norm_2( p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(),
                              p_element->GetNode(1)->rGetLocation()) );

        linear_force.SetCutOffLength(dist-0.1);

        // Coverage
        TS_ASSERT_DELTA(linear_force.GetCutOffLength(), dist-0.1, 1e-4);

        force_on_spring = linear_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),
                                                                  p_element->GetNodeGlobalIndex(0),
                                                                  cell_population);
        TS_ASSERT_DELTA(force_on_spring[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);
    }

    void TestGeneralisedLinearSpringForceCalculationIn1d() throw (Exception)
    {
        // Create a 1D mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes(), std::vector<unsigned>(), DIFFERENTIATED);

        // Create cell population
        std::vector<CellPtr> cells_copy(cells);
        MeshBasedCellPopulation<1> cell_population(mesh, cells);

        // Create force law object
        GeneralisedLinearSpringForce<1> linear_force;

        // Initialise a vector of node forces
        std::vector<c_vector<double, 1> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(1));
        }

        // Compute forces on nodes
        linear_force.AddForceContribution(node_forces, cell_population);

        // Test that all springs are in equilibrium
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(node_forces[node_index](0), 0.0, 1e-6);
        }

        // Scale entire mesh and check that forces are correctly calculated
        double scale_factor = 1.5;
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            c_vector<double,1> old_point = mesh.GetNode(node_index)->rGetLocation();
            ChastePoint<1> new_point;
            new_point.rGetLocation()[0] = scale_factor*old_point[0];
            mesh.SetNode(node_index, new_point, false);
        }

        // Recalculate node forces (we can re-use node_forces
        // as previously each node had zero net force on it)
        linear_force.AddForceContribution(node_forces, cell_population);

        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            if (node_index == 0)
            {
                // The first node only experiences a force from its neighbour to the right
                TS_ASSERT_DELTA(node_forces[node_index](0), linear_force.GetMeinekeSpringStiffness()*(scale_factor-1), 1e-6);
            }
            else if (node_index == cell_population.GetNumNodes()-1)
            {
                // The last node only experiences a force from its neighbour to the left
                TS_ASSERT_DELTA(node_forces[node_index](0), -linear_force.GetMeinekeSpringStiffness()*(scale_factor-1), 1e-6);
            }
            else
            {
                // The net force on each interior node should still be zero
                TS_ASSERT_DELTA(node_forces[node_index](0), 0.0, 1e-6);
            }
        }

        // Create another cell population and force law
        MutableMesh<1,1> mesh2;
        mesh2.ConstructLinearMesh(5);

        MeshBasedCellPopulation<1> cell_population2(mesh2, cells_copy);
        GeneralisedLinearSpringForce<1> linear_force2;

        // Move one node and check that forces are correctly calculated
        ChastePoint<1> shifted_point;
        shifted_point.rGetLocation()[0] = 2.5;
        mesh2.SetNode(2, shifted_point);

        c_vector<double,1> force_between_1_and_2 = linear_force2.CalculateForceBetweenNodes(1, 2, cell_population2);
        TS_ASSERT_DELTA(force_between_1_and_2[0], linear_force.GetMeinekeSpringStiffness()*0.5, 1e-6);

        c_vector<double,1> force_between_2_and_3 = linear_force2.CalculateForceBetweenNodes(2, 3, cell_population2);
        TS_ASSERT_DELTA(force_between_2_and_3[0], -linear_force.GetMeinekeSpringStiffness()*0.5, 1e-6);

        // Initialise a vector of node forces
        std::vector<c_vector<double,1> > node_forces2;
        node_forces2.reserve(cell_population2.GetNumNodes());

        for (unsigned i=0; i<cell_population2.GetNumNodes(); i++)
        {
             node_forces2.push_back(zero_vector<double>(1));
        }

        linear_force2.AddForceContribution(node_forces2, cell_population2);

        TS_ASSERT_DELTA(node_forces2[2](0), -linear_force.GetMeinekeSpringStiffness(), 1e-6);
    }

    void TestGeneralisedLinearSpringForceCalculationIn3d() throw (Exception)
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-50.0);
        }

        std::vector<CellPtr> cells_copy(cells);
        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        GeneralisedLinearSpringForce<3> linear_force;

        // Test forces on springs
        unsigned nodeA = 0, nodeB = 1;
        Element<3,3>* p_element = mesh.GetElement(0);
        c_vector<double, 3> force = linear_force.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(nodeA),
                                                                            p_element->GetNodeGlobalIndex(nodeB),
                                                                            cell_population);
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(force[i], 0.0, 1e-6);
        }

        // Initialise a vector of node forces
        std::vector<c_vector<double, 3> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(3));
        }

        linear_force.AddForceContribution(node_forces, cell_population);

        for (unsigned j=0; j<4; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(node_forces[j](k), 0.0, 1e-6);
            }
        }

        // Scale entire mesh and check that forces are correctly calculated
        double scale_factor = 1.5;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> old_point = mesh.GetNode(i)->rGetLocation();
            ChastePoint<3> new_point;
            new_point.rGetLocation()[0] = scale_factor*old_point[0];
            new_point.rGetLocation()[1] = scale_factor*old_point[1];
            new_point.rGetLocation()[2] = scale_factor*old_point[2];
            mesh.SetNode(i, new_point, false);
        }

        // Recalculate node forces (we can just re-use node_forces,
        // as previously each node had zero net force on it)
        linear_force.AddForceContribution(node_forces, cell_population);

        for (unsigned j=0; j<4; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(fabs(node_forces[j](k)), linear_force.GetMeinekeSpringStiffness()*(scale_factor-1)*sqrt(2),1e-6);
            }
        }

        // Move one node and check that forces are correctly calculated
        MutableMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);

        MeshBasedCellPopulation<3> cell_population2(mesh2, cells_copy);
        GeneralisedLinearSpringForce<3> linear_force2;

        c_vector<double,3> old_point = mesh2.GetNode(0)->rGetLocation();
        ChastePoint<3> new_point;
        new_point.rGetLocation()[0] = 0.0;
        new_point.rGetLocation()[1] = 0.0;
        new_point.rGetLocation()[2] = 0.0;
        mesh2.SetNode(0, new_point, false);

        unsigned nodeA2 = 0, nodeB2 = 1;
        Element<3,3>* p_element2 = mesh2.GetElement(0);
        c_vector<double,3> force2 = linear_force2.CalculateForceBetweenNodes(p_element2->GetNodeGlobalIndex(nodeA2),
                                                                             p_element2->GetNodeGlobalIndex(nodeB2),
                                                                             cell_population2);

        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(fabs(force2[i]),linear_force.GetMeinekeSpringStiffness()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }

        // Initialise a vector of node forces
        std::vector<c_vector<double,3> > node_forces2;
        node_forces2.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             node_forces2.push_back(zero_vector<double>(3));
        }

        linear_force2.AddForceContribution(node_forces2, cell_population2);

        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(node_forces2[0](i),linear_force.GetMeinekeSpringStiffness()*(1 - sqrt(3)/(2*sqrt(2)))/sqrt(3.0),1e-6);
        }
    }

    void TestForceOutputParameters()
    {
        std::string output_directory = "TestForcesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with GeneralisedLinearSpringForce
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(linear_force.GetIdentifier(), "GeneralisedLinearSpringForce-2");

        out_stream linear_force_parameter_file = output_file_handler.OpenOutputFile("linear_results.parameters");
        linear_force.OutputForceParameters(linear_force_parameter_file);
        linear_force_parameter_file->close();

        std::string linear_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + linear_force_results_dir + "linear_results.parameters cell_based/test/data/TestForces/linear_results.parameters").c_str()), 0);

        // Test with ChemotacticForce
        ChemotacticForce<2> chemotactic_force;
        TS_ASSERT_EQUALS(chemotactic_force.GetIdentifier(), "ChemotacticForce-2");

        out_stream chemotactic_force_parameter_file = output_file_handler.OpenOutputFile("chemotactic_results.parameters");
        chemotactic_force.OutputForceParameters(chemotactic_force_parameter_file);
        chemotactic_force_parameter_file->close();

        std::string chemotactic_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + chemotactic_force_results_dir + "chemotactic_results.parameters cell_based/test/data/TestForces/chemotactic_results.parameters").c_str()), 0);

        // Test with RepulsionForce
        RepulsionForce<2> repulsion_force;
        TS_ASSERT_EQUALS(repulsion_force.GetIdentifier(), "RepulsionForce-2");

        out_stream repulsion_force_parameter_file = output_file_handler.OpenOutputFile("repulsion_results.parameters");
        repulsion_force.OutputForceParameters(repulsion_force_parameter_file);
        repulsion_force_parameter_file->close();

        std::string repulsion_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + repulsion_force_results_dir + "repulsion_results.parameters cell_based/test/data/TestForces/repulsion_results.parameters").c_str()), 0);

        // Test with NagaiHondaForce
        NagaiHondaForce<2> nagai_force;
        TS_ASSERT_EQUALS(nagai_force.GetIdentifier(), "NagaiHondaForce-2");

        out_stream nagai_force_parameter_file = output_file_handler.OpenOutputFile("nagai_results.parameters");
        nagai_force.OutputForceParameters(nagai_force_parameter_file);
        nagai_force_parameter_file->close();

        std::string nagai_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + nagai_force_results_dir + "nagai_results.parameters cell_based/test/data/TestForces/nagai_results.parameters").c_str()), 0);

        // Test with WelikyOsterForce
        WelikyOsterForce<2> weliky_force;
        TS_ASSERT_EQUALS(weliky_force.GetIdentifier(), "WelikyOsterForce-2");

        out_stream weliky_force_parameter_file = output_file_handler.OpenOutputFile("weliky_results.parameters");
        weliky_force.OutputForceParameters(weliky_force_parameter_file);
        weliky_force_parameter_file->close();

        std::string weliky_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + weliky_force_results_dir + "weliky_results.parameters cell_based/test/data/TestForces/weliky_results.parameters").c_str()), 0);

        // Test with DiffusionForce
		DiffusionForce<2> diffusion_force;
		TS_ASSERT_EQUALS(diffusion_force.GetIdentifier(), "DiffusionForce-2");

		out_stream diffusion_force_parameter_file = output_file_handler.OpenOutputFile("diffusion_results.parameters");
		diffusion_force.OutputForceParameters(diffusion_force_parameter_file);
		diffusion_force_parameter_file->close();

		std::string diffusion_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
		TS_ASSERT_EQUALS(system(("diff " + diffusion_force_results_dir + "diffusion_results.parameters cell_based/test/data/TestForces/diffusion_results.parameters").c_str()), 0);
    }

    void TestGeneralisedLinearSpringForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "GeneralisedLinearSpringForce.arch";

        {
            GeneralisedLinearSpringForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetMeinekeSpringStiffness(12.34);
            force.SetMeinekeDivisionRestingSpringLength(0.856);
            force.SetMeinekeSpringGrowthDuration(2.593);

            // Serialize via pointer to most abstract class possible
            AbstractForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;

            // Test member variables
            TS_ASSERT_DELTA((static_cast<GeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringStiffness(), 12.34, 1e-6);
            TS_ASSERT_DELTA((static_cast<GeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
            TS_ASSERT_DELTA((static_cast<GeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestChemotacticForceMethods() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 5;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MAKE_PTR(CellLabel, p_label);
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-10);
            cells[i]->AddCellProperty(p_label);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cellwise data and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetPopulationAndNumVars(&cell_population, 1);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            p_data->SetValue(x/50.0, p_mesh->GetNode(i)->GetIndex());
        }

        ChemotacticForce<2> chemotactic_force;

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             node_forces.push_back(zero_vector<double>(2));
        }
        chemotactic_force.AddForceContribution(node_forces, cell_population);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
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
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ChemotacticForce.arch";

        {
            ChemotacticForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // No member variables to set

            // Serialize via pointer to most abstract class possible
            AbstractForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;

            // No member variables to test, so just test a method
            TS_ASSERT_DELTA((static_cast<ChemotacticForce<2>*>(p_force))->GetChemotacticForceMagnitude(12.0, 3.5), 12.0, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestRepulsionForceMethods() throw (Exception)
    {
        // Create a NodeBasedCellPopulation
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(2, true, 3.0, 0.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(100.0);
        cell_population.Update(); //Needs to be called seperately as not in a simulation

        RepulsionForce<2> repulsion_force;

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }
        repulsion_force.AddForceContribution(node_forces, cell_population);

        // First two cells repel each other and second 2 cells are too far appart
        TS_ASSERT_DELTA(node_forces[0][0], -34.5387, 1e-4);
        TS_ASSERT_DELTA(node_forces[0][1], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_forces[1][0], 34.5387, 1e-4);
        TS_ASSERT_DELTA(node_forces[1][1], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_forces[2][0], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_forces[2][1], 0.0, 1e-4);

        // When the mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestRepulsionForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "RepulsionForce.arch";

        {
            RepulsionForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // No extra member variables, so set member variables on parent class
            force.SetMeinekeSpringStiffness(12.35);
            force.SetMeinekeDivisionRestingSpringLength(0.756);
            force.SetMeinekeSpringGrowthDuration(2.693);

            // Serialize via pointer to most abstract class possible
            AbstractForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;

            // No extra member variables, so test member variables on parent class
            TS_ASSERT_DELTA((static_cast<RepulsionForce<2>*>(p_force))->GetMeinekeSpringStiffness(), 12.35, 1e-6);
            TS_ASSERT_DELTA((static_cast<RepulsionForce<2>*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.756, 1e-6);
            TS_ASSERT_DELTA((static_cast<RepulsionForce<2>*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.693, 1e-6);

            // Tidy up
            delete p_force;
        }
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
            nodes.push_back(new Node<2>(i, true, cos(angles[i]), sin(angles[i])));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Set up the cell
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(DIFFERENTIATED);

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetBirthTime(-1.0);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.InitialiseCells();

        // Create a force system
        NagaiHondaForce<2> force;

        // Test get/set methods
        TS_ASSERT_DELTA(force.GetNagaiHondaDeformationEnergyParameter(), 100.0, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaMembraneSurfaceEnergyParameter(), 10.0, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaCellCellAdhesionEnergyParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaCellBoundaryAdhesionEnergyParameter(), 1.0, 1e-12);

        force.SetNagaiHondaDeformationEnergyParameter(5.8);
        force.SetNagaiHondaMembraneSurfaceEnergyParameter(17.9);
        force.SetNagaiHondaCellCellAdhesionEnergyParameter(0.5);
        force.SetNagaiHondaCellBoundaryAdhesionEnergyParameter(0.6);

        TS_ASSERT_DELTA(force.GetNagaiHondaDeformationEnergyParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaMembraneSurfaceEnergyParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaCellCellAdhesionEnergyParameter(), 0.5, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-12);

        force.SetNagaiHondaDeformationEnergyParameter(100.0);
        force.SetNagaiHondaMembraneSurfaceEnergyParameter(10.0);
        force.SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        force.SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        force.AddForceContribution(node_forces, cell_population);

        // The force on each node should be radially inward, with the same magnitude for all nodes
        double force_magnitude = norm_2(node_forces[0]);

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(node_forces[i]), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(node_forces[i][0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(node_forces[i][1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        double normal_target_area = force.GetTargetAreaOfCell(cell_population.GetCellUsingLocationIndex(0));

        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.25, 2);

        // Set the cell to be necrotic
        cell_population.GetCellUsingLocationIndex(0)->StartApoptosis();

        double initial_apoptotic_target_area = force.GetTargetAreaOfCell(cell_population.GetCellUsingLocationIndex(0));

        TS_ASSERT_DELTA(normal_target_area, initial_apoptotic_target_area, 1e-6);

        // Reset force vector
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces[i] = zero_vector<double>(2);
        }

        force.AddForceContribution(node_forces, cell_population);

        // The force on each node should not yet be affected by setting the cell to be apoptotic
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(node_forces[i]), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(node_forces[i][0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(node_forces[i][1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Increment time
        p_simulation_time->IncrementTimeOneStep();

        double later_apoptotic_target_area = force.GetTargetAreaOfCell(cell_population.GetCellUsingLocationIndex(0));

        TS_ASSERT_LESS_THAN(later_apoptotic_target_area, initial_apoptotic_target_area);

        TS_ASSERT_DELTA(cell_population.GetCellUsingLocationIndex(0)->GetTimeUntilDeath(), 0.125, 1e-6);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces[i] = zero_vector<double>(2);
        }

        force.AddForceContribution(node_forces, cell_population);

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

    void TestNagaiHondaGetTargetAreaOfCell() throw (Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3*0.25, 3);

        // Create mesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        cells[0]->GetCellCycleModel()->SetCellProliferativeType(DIFFERENTIATED);
        cells[4]->GetCellCycleModel()->SetCellProliferativeType(DIFFERENTIATED);
        for (unsigned i=0; i<cells.size(); i++)
        {
            double birth_time = 0.0 - 2*i;
            cells[i]->SetBirthTime(birth_time);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.InitialiseCells(); // this method must be called explicitly as there is no simulation

        // Create a force system
        NagaiHondaForce<2> force;

        // Check GetTargetAreaOfCell()
        for (VertexMesh<2,2>::VertexElementIterator iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(elem_index);
            double expected_area = force.GetMatureCellTargetArea();

            if (elem_index!=4 && elem_index<=7)
            {
                expected_area *= 0.5*(1 + ((double)elem_index)/7.0);
            }

            double actual_area = force.GetTargetAreaOfCell(p_cell);

            TS_ASSERT_DELTA(actual_area, expected_area, 1e-12);
        }

        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        CellPtr p_cell_1 = cell_population.GetCellUsingLocationIndex(1);
        CellPtr p_cell_4 = cell_population.GetCellUsingLocationIndex(4);

        // Make cell 1 and 4 undergo apoptosis
        p_cell_1->StartApoptosis();
        p_cell_4->StartApoptosis();

        double actual_area_0 = force.GetTargetAreaOfCell(p_cell_0);
        double actual_area_1 = force.GetTargetAreaOfCell(p_cell_1);
        double actual_area_4 = force.GetTargetAreaOfCell(p_cell_4);

        double expected_area_0 = 0.5;
        double expected_area_1 = force.GetMatureCellTargetArea()*0.5*(1.0 + 1.0/7.0);
        double expected_area_4 = force.GetMatureCellTargetArea();

        TS_ASSERT_DELTA(actual_area_0, expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1, expected_area_1, 1e-12);
        TS_ASSERT_DELTA(actual_area_4, expected_area_4, 1e-12);

        // Run for one time step
        p_simulation_time->IncrementTimeOneStep();

        double actual_area_0_after_dt = force.GetTargetAreaOfCell(p_cell_0);
        double actual_area_1_after_dt = force.GetTargetAreaOfCell(p_cell_1);
        double actual_area_4_after_dt = force.GetTargetAreaOfCell(p_cell_4);

        // The target areas of cells 1 and 4 should have halved
        expected_area_0 = force.GetMatureCellTargetArea()*0.5*(1.0 + 0.5*0.25);

        TS_ASSERT_DELTA(actual_area_0_after_dt, expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_dt, 0.5*expected_area_1, 1e-12);
        TS_ASSERT_DELTA(actual_area_4_after_dt, 0.5*expected_area_4, 1e-12);

        // Make cell 0 undergo apoptosis
        p_cell_0->StartApoptosis();

        // Now run on for a further time step
        p_simulation_time->IncrementTimeOneStep();

        double actual_area_0_after_2dt = force.GetTargetAreaOfCell(p_cell_0);
        double actual_area_1_after_2dt = force.GetTargetAreaOfCell(p_cell_1);
        double actual_area_4_after_2dt = force.GetTargetAreaOfCell(p_cell_4);

        // Cells 1 and 4 should now have zero target area and the target area of cell 0 should have halved
        TS_ASSERT_DELTA(actual_area_0_after_2dt, 0.5*expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_2dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_4_after_2dt, 0.0, 1e-12);

        // Now run on for even further, for coverage
        p_simulation_time->IncrementTimeOneStep();

        double actual_area_0_after_3dt = force.GetTargetAreaOfCell(p_cell_0);
        double actual_area_1_after_3dt = force.GetTargetAreaOfCell(p_cell_1);
        double actual_area_4_after_3dt = force.GetTargetAreaOfCell(p_cell_4);

        // All apoptotic cells should now have zero target area
        TS_ASSERT_DELTA(actual_area_0_after_3dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_3dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_4_after_3dt, 0.0, 1e-12);
    }

    void TestNagaiHondaForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "NagaiHondaForce.arch";

        {
            NagaiHondaForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetNagaiHondaDeformationEnergyParameter(5.8);
            force.SetNagaiHondaMembraneSurfaceEnergyParameter(17.9);
            force.SetNagaiHondaCellCellAdhesionEnergyParameter(0.5);
            force.SetNagaiHondaCellBoundaryAdhesionEnergyParameter(0.6);
            force.SetMatureCellTargetArea(0.7);

            // Serialize via pointer to most abstract class possible
            AbstractForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(static_cast<NagaiHondaForce<2>*>(p_force)->GetNagaiHondaDeformationEnergyParameter(), 5.8, 1e-12);
            TS_ASSERT_DELTA(static_cast<NagaiHondaForce<2>*>(p_force)->GetNagaiHondaMembraneSurfaceEnergyParameter(), 17.9, 1e-12);
            TS_ASSERT_DELTA(static_cast<NagaiHondaForce<2>*>(p_force)->GetNagaiHondaCellCellAdhesionEnergyParameter(), 0.5, 1e-12);
            TS_ASSERT_DELTA(static_cast<NagaiHondaForce<2>*>(p_force)->GetNagaiHondaCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-12);
            TS_ASSERT_DELTA(static_cast<NagaiHondaForce<2>*>(p_force)->GetMatureCellTargetArea(), 0.7, 1e-12);

            // Tidy up
            delete p_force;
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
            nodes.push_back(new Node<2>(i, true, cos(angles[i]), sin(angles[i])));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Set up the cell
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(DIFFERENTIATED);

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetBirthTime(-1.0);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.InitialiseCells();

        // Create a force system
        WelikyOsterForce<2> force;

        // Test set/get methods
        TS_ASSERT_DELTA(force.GetWelikyOsterAreaParameter(), 1.0, 1e-6);
        TS_ASSERT_DELTA(force.GetWelikyOsterPerimeterParameter(), 1.0, 1e-6);

        force.SetWelikyOsterAreaParameter(15.0);
        force.SetWelikyOsterPerimeterParameter(17.0);

        TS_ASSERT_DELTA(force.GetWelikyOsterAreaParameter(), 15.0, 1e-6);
        TS_ASSERT_DELTA(force.GetWelikyOsterPerimeterParameter(), 17.0, 1e-6);

        force.SetWelikyOsterAreaParameter(1.0);
        force.SetWelikyOsterPerimeterParameter(1.0);

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        force.AddForceContribution(node_forces, cell_population);

        // The force on each node should be radially inward, with the same magnitude for all nodes
        double force_magnitude = norm_2(node_forces[0]);

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(node_forces[i]), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(node_forces[i][0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(node_forces[i][1], -force_magnitude*sin(angles[i]), 1e-4);
        }
    }

    void TestWelikyOsterForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "WelikyOsterForce.arch";

        {
            WelikyOsterForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetWelikyOsterAreaParameter(15.12);
            force.SetWelikyOsterPerimeterParameter(17.89);

            // Serialize via pointer to most abstract class possible
            AbstractForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(static_cast<WelikyOsterForce<2>*>(p_force)->GetWelikyOsterAreaParameter(), 15.12, 1e-12);
            TS_ASSERT_DELTA(static_cast<WelikyOsterForce<2>*>(p_force)->GetWelikyOsterPerimeterParameter(), 17.89, 1e-12);

            // Tidy up
            delete p_force;
        }
    }

    void TestCentreBasedForcesWithVertexCellPopulation() throw (Exception)
    {
        // Construct simple vertex mesh
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

        MutableVertexMesh<2,2> mesh(nodes, elements, 0.01, 2.0);

        // Create cell
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetBirthTime(-1.0);
        cells.push_back(p_cell);

        // Create VertexBasedCellPopulation
        VertexBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.InitialiseCells();

        // Create a vector forces on nodes
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(num_nodes);
        for (unsigned i=0; i<num_nodes; i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        // Test that a subclass of AbstractTwoBodyInteractionForce throws the correct exception
        GeneralisedLinearSpringForce<2> spring_force;
        TS_ASSERT_THROWS_THIS(spring_force.AddForceContribution(node_forces, cell_population),
                "Subclasses of AbstractTwoBodyInteractionForce are to be used with subclasses of AbstractCentreBasedCellPopulation only");

        // Test that RepulsionForce throws the correct exception
        RepulsionForce<2> repulsion_force;
        TS_ASSERT_THROWS_THIS(repulsion_force.AddForceContribution(node_forces, cell_population),
                 "RepulsionForce is to be used with a NodeBasedCellPopulation only");
    }

    void TestIncorrectForcesWithNodeBasedCellPopulation() throw (Exception)
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

        // Test that NagaiHondaForce throws the correct exception
        NagaiHondaForce<2> nagai_honda_force;
        TS_ASSERT_THROWS_THIS(nagai_honda_force.AddForceContribution(node_forces, cell_population),
                "NagaiHondaForce is to be used with a VertexBasedCellPopulation only");

        // Test that WelikyOsterForce throws the correct exception
        WelikyOsterForce<2> weliky_oster_force;
        TS_ASSERT_THROWS_THIS(weliky_oster_force.AddForceContribution(node_forces, cell_population),
                "WelikyOsterForce is to be used with a VertexBasedCellPopulation only");

        // When the node-only mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestDiffusionForceIn1D()
    {
        // Set up time parameters
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create a 1D mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes(), std::vector<unsigned>(), DIFFERENTIATED);

        // Create cell population
        std::vector<CellPtr> cells_copy(cells);
        MeshBasedCellPopulation<1> cell_population(mesh, cells);

        // Create force law object
        DiffusionForce<1> diffusion_force;

        // Initialise a vector of node forces
        std::vector<c_vector<double, 1> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(1));
        }

        // Compute forces on nodes
        diffusion_force.AddForceContribution(node_forces, cell_population);

        // Test Set and Get methods for the diffusion force
        TS_ASSERT_DELTA(diffusion_force.GetCutOffLength(), 10.0, 1e-10);
        TS_ASSERT_DELTA(diffusion_force.GetDiffusionConstant(), 0.01, 1e-10);
        diffusion_force.SetCutOffLength(20.0);
        diffusion_force.SetDiffusionConstant(0.1);
        TS_ASSERT_DELTA(diffusion_force.GetCutOffLength(), 20.0, 1e-10);
        TS_ASSERT_DELTA(diffusion_force.GetDiffusionConstant(), 0.1, 1e-10);
        diffusion_force.SetCutOffLength(10.0);
        diffusion_force.SetDiffusionConstant(0.01);
    }

    void TestDiffusionForceIn2D()
    {
        // Define the seed
        RandomNumberGenerator::Instance()->Reseed(0);

        // Set up time parameters
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create a NodeBasedCellPopulation
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(100.0);
        cell_population.Update(); //Needs to be called separately as not in a simulation

        // Create force law object
        DiffusionForce<2> diffusion_force;

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        double variance=0.0;
        // Loop over time iterations
        unsigned number_iteration=1000u;
        for (unsigned i=0; i<number_iteration; i++)
        {
            // Re-initialize the force on node zero
            node_forces[0]=zero_vector<double>(2);

            // Compute forces on nodes
            diffusion_force.AddForceContribution(node_forces, cell_population);

            // Calculate the variance
            variance+=pow(norm_2(node_forces[0]),2);
        }

        int dim = 2;
        variance /= number_iteration*2*dim*diffusion_force.GetDiffusionConstant()*SimulationTime::Instance()->GetTimeStep();
        TS_ASSERT_DELTA(variance, 1.0, 1e-1)

        // When the node-only mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

        // Tidy up
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    void TestDiffusionForceIn3D()
    {
        // Define the seed
        RandomNumberGenerator::Instance()->Reseed(0);

        // Set up time parameters
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create a NodeBasedCellPopulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(100.0);
        cell_population.Update(); //Needs to be called separately as not in a simulation

        // Create force law object
        DiffusionForce<3> diffusion_force;

        // Initialise a vector of node forces
        std::vector<c_vector<double, 3> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(3));
        }

        double variance=0.0;
        // Loop over time iterations
        unsigned number_iteration=1000u;
        for (unsigned i=0; i<number_iteration; i++)
        {
            // Re-initialize the force on node zero
            node_forces[0]=zero_vector<double>(3);

            // Compute forces on nodes
            diffusion_force.AddForceContribution(node_forces, cell_population);

            // Calculate the variance
            variance+=pow(norm_2(node_forces[0]),2);
        }

        int dim=3;
        variance /= number_iteration*2*dim*diffusion_force.GetDiffusionConstant()*SimulationTime::Instance()->GetTimeStep();
        TS_ASSERT_DELTA(variance, 1.0, 1e-1)

        // When the node-only mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

        // Tidy up
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    void TestDiffusionForceArchiving() throw (Exception)
        {
            OutputFileHandler handler("archive", false);
            std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DiffusionForce.arch";

            {
                DiffusionForce<2> force;

                std::ofstream ofs(archive_filename.c_str());
                boost::archive::text_oarchive output_arch(ofs);

                // Serialize via pointer to most abstract class possible
                AbstractForce<2>* const p_force = &force;
                output_arch << p_force;
            }

            {
                AbstractForce<2>* p_force;

                // Create an input archive
                std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
                boost::archive::text_iarchive input_arch(ifs);

                // Restore from the archive
                input_arch >> p_force;

                // Test member variables
                TS_ASSERT_DELTA((static_cast<DiffusionForce<2>*>(p_force))->GetCutOffLength(), 10.0, 1e-6);
                TS_ASSERT_DELTA((static_cast<DiffusionForce<2>*>(p_force))->GetDiffusionConstant(), 0.01, 1e-6);

                // Tidy up
                delete p_force;
            }
        }
};

#endif /*TESTFORCES_HPP_*/
