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

#ifndef TESTFORCES_HPP_
#define TESTFORCES_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ChemotacticForce.hpp"
#include "RepulsionForce.hpp"
#include "NagaiHondaForce.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "WelikyOsterForce.hpp"
#include "FarhadifarForce.hpp"
#include "DiffusionForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestForces : public AbstractCellBasedTestSuite
{
public:

    void TestGeneralisedLinearSpringForceMethods()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
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
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Test node force calculation
        linear_force.AddForceContribution(cell_population);

        // Test forces on non-ghost nodes
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
        }

        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point;
        old_point = p_mesh->GetNode(59)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];

        p_mesh->SetNode(59, new_point, false);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        linear_force.AddForceContribution(cell_population);

        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);

        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], (-3+4.0/sqrt(7.0))*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);

        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 0.5*linear_force.GetMeinekeSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);

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

    void TestGeneralisedLinearSpringForceCalculationIn1d()
    {
        // Create a 1D mesh with nodes equally spaced a unit distance apart
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes(), std::vector<unsigned>(), p_diff_type);

        // Create cell population
        std::vector<CellPtr> cells_copy(cells);
        MeshBasedCellPopulation<1> cell_population(mesh, cells);

        // Create force law object
        GeneralisedLinearSpringForce<1> linear_force;

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Compute forces on nodes
        linear_force.AddForceContribution(cell_population);

        // Test that all springs are in equilibrium
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0, 1e-6);
        }

        // Scale entire mesh and check that forces are correctly calculated
        double scale_factor = 1.5;
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
            c_vector<double,1> old_point;
            old_point = mesh.GetNode(node_index)->rGetLocation();

            ChastePoint<1> new_point;
            new_point.rGetLocation()[0] = scale_factor*old_point[0];
            mesh.SetNode(node_index, new_point, false);
        }

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        linear_force.AddForceContribution(cell_population);

        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            if (node_index == 0)
            {
                // The first node only experiences a force from its neighbour to the right
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], linear_force.GetMeinekeSpringStiffness()*(scale_factor-1), 1e-6);
            }
            else if (node_index == cell_population.GetNumNodes()-1)
            {
                // The last node only experiences a force from its neighbour to the left
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], -linear_force.GetMeinekeSpringStiffness()*(scale_factor-1), 1e-6);
            }
            else
            {
                // The net force on each interior node should still be zero
                TS_ASSERT_DELTA(cell_population.GetNode(node_index)->rGetAppliedForce()[0], 0.0, 1e-6);
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

        for (unsigned i=0; i<cell_population2.GetNumNodes(); i++)
        {
             cell_population2.GetNode(i)->ClearAppliedForce();
        }

        linear_force2.AddForceContribution(cell_population2);

        TS_ASSERT_DELTA(cell_population2.GetNode(2)->rGetAppliedForce()[0], -linear_force.GetMeinekeSpringStiffness(), 1e-6);
    }

    void TestGeneralisedLinearSpringForceCalculationIn3d()
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
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

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        linear_force.AddForceContribution(cell_population);

        for (unsigned j=0; j<4; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(cell_population.GetNode(j)->rGetAppliedForce()[k], 0.0, 1e-6);
            }
        }

        // Scale entire mesh and check that forces are correctly calculated
        double scale_factor = 1.5;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> old_point;
            old_point = mesh.GetNode(i)->rGetLocation();
            ChastePoint<3> new_point;
            new_point.rGetLocation()[0] = scale_factor*old_point[0];
            new_point.rGetLocation()[1] = scale_factor*old_point[1];
            new_point.rGetLocation()[2] = scale_factor*old_point[2];
            mesh.SetNode(i, new_point, false);
        }

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        linear_force.AddForceContribution(cell_population);

        for (unsigned j=0; j<4; j++)
        {
            for (unsigned k=0; k<3; k++)
            {
                TS_ASSERT_DELTA(fabs(cell_population.GetNode(j)->rGetAppliedForce()[k]), linear_force.GetMeinekeSpringStiffness()*(scale_factor-1)*sqrt(2.0),1e-6);
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
            TS_ASSERT_DELTA(fabs(force2[i]),linear_force.GetMeinekeSpringStiffness()*(1 - sqrt(3.0)/(2*sqrt(2.0)))/sqrt(3.0),1e-6);
        }

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population2.GetNode(i)->ClearAppliedForce();
        }

        linear_force2.AddForceContribution(cell_population2);

        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(cell_population2.GetNode(0)->rGetAppliedForce()[i],linear_force.GetMeinekeSpringStiffness()*(1 - sqrt(3.0)/(2*sqrt(2.0)))/sqrt(3.0),1e-6);
        }
    }

    void TestDifferentialAdhesionGeneralisedLinearSpringForceMethods()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size(), location_indices);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create force
        DifferentialAdhesionGeneralisedLinearSpringForce<2> force;

        // Test set/get method
        TS_ASSERT_DELTA(force.GetHomotypicLabelledSpringConstantMultiplier(), 1.0, 1e-6);
        TS_ASSERT_DELTA(force.GetHeterotypicSpringConstantMultiplier(), 1.0, 1e-6);

        force.SetHomotypicLabelledSpringConstantMultiplier(2.0);
        force.SetHeterotypicSpringConstantMultiplier(4.0);

        TS_ASSERT_DELTA(force.GetHomotypicLabelledSpringConstantMultiplier(), 2.0, 1e-6);
        TS_ASSERT_DELTA(force.GetHeterotypicSpringConstantMultiplier(), 4.0, 1e-6);

        // Initialise a vector of node forces
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point;
        old_point = p_mesh->GetNode(59)->rGetLocation();
        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];
        p_mesh->SetNode(59, new_point, false);

        double spring_stiffness = force.GetMeinekeSpringStiffness();

        // Test the case where node 59 and its neighbours are unlabelled
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        force.AddForceContribution(cell_population);

        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], (-3+4.0/sqrt(7.0))*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);

        // Next, test the case where node 59 is labelled but its neighbours are not...
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
        cell_population.GetCellUsingLocationIndex(59)->AddCellProperty(p_label);

        force.AddForceContribution(cell_population);

        // ...for which the force magnitude should be increased by 4, our chosen multiplier for heterotypic interactions under attraction
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 4.0*0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], -0.5*spring_stiffness+4.0*(4.0/sqrt(7.0)-2.5)*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);

        // Finally, test the case where node 59 and its neighbours are labelled...
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->AddCellProperty(p_label);
        }

        force.AddForceContribution(cell_population);

        // ...for which the force magnitude should be increased by 2, our chosen multiplier for homotypic labelled interactions, again only for attractive interactions
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 2.0*0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], -0.5*spring_stiffness+2.0*(4.0/sqrt(7.0)-2.5)*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);
    }

    void TestForceOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestForcesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with GeneralisedLinearSpringForce
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(linear_force.GetIdentifier(), "GeneralisedLinearSpringForce-2-2");

        out_stream linear_force_parameter_file = output_file_handler.OpenOutputFile("linear_results.parameters");
        linear_force.OutputForceParameters(linear_force_parameter_file);
        linear_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("linear_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/linear_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with DifferentialAdhesionGeneralisedLinearSpringForce
        DifferentialAdhesionGeneralisedLinearSpringForce<2> differential_linear_force;
        differential_linear_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(differential_linear_force.GetIdentifier(), "DifferentialAdhesionGeneralisedLinearSpringForce-2-2");

        out_stream differential_linear_force_parameter_file = output_file_handler.OpenOutputFile("differential_linear_results.parameters");
        differential_linear_force.OutputForceParameters(differential_linear_force_parameter_file);
        differential_linear_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("differential_linear_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/differential_linear_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with ChemotacticForce
        ChemotacticForce<2> chemotactic_force;
        TS_ASSERT_EQUALS(chemotactic_force.GetIdentifier(), "ChemotacticForce-2");

        out_stream chemotactic_force_parameter_file = output_file_handler.OpenOutputFile("chemotactic_results.parameters");
        chemotactic_force.OutputForceParameters(chemotactic_force_parameter_file);
        chemotactic_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("chemotactic_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/chemotactic_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with RepulsionForce
        RepulsionForce<2> repulsion_force;
        TS_ASSERT_EQUALS(repulsion_force.GetIdentifier(), "RepulsionForce-2");

        out_stream repulsion_force_parameter_file = output_file_handler.OpenOutputFile("repulsion_results.parameters");
        repulsion_force.OutputForceParameters(repulsion_force_parameter_file);
        repulsion_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("repulsion_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/repulsion_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with NagaiHondaForce
        NagaiHondaForce<2> nagai_force;
        TS_ASSERT_EQUALS(nagai_force.GetIdentifier(), "NagaiHondaForce-2");

        out_stream nagai_force_parameter_file = output_file_handler.OpenOutputFile("nagai_results.parameters");
        nagai_force.OutputForceParameters(nagai_force_parameter_file);
        nagai_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("nagai_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/nagai_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with NagaiHondaDifferentialAdhesionForce
        NagaiHondaDifferentialAdhesionForce<2> nagai_da_force;
        TS_ASSERT_EQUALS(nagai_da_force.GetIdentifier(), "NagaiHondaDifferentialAdhesionForce-2");

        out_stream nagai_da_force_parameter_file = output_file_handler.OpenOutputFile("nagai_da_results.parameters");
        nagai_da_force.OutputForceParameters(nagai_force_parameter_file);
        nagai_da_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("nagai_da_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/nagai_da_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with WelikyOsterForce
        WelikyOsterForce<2> weliky_force;
        TS_ASSERT_EQUALS(weliky_force.GetIdentifier(), "WelikyOsterForce-2");

        out_stream weliky_force_parameter_file = output_file_handler.OpenOutputFile("weliky_results.parameters");
        weliky_force.OutputForceParameters(weliky_force_parameter_file);
        weliky_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("weliky_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/weliky_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        FarhadifarForce<2> force;
        TS_ASSERT_EQUALS(force.GetIdentifier(), "FarhadifarForce-2");

        out_stream farhadifar_force_parameter_file = output_file_handler.OpenOutputFile("farhadifar_results.parameters");
        force.OutputForceParameters(farhadifar_force_parameter_file);
        farhadifar_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("farhadifar_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/farhadifar_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with DiffusionForce
        DiffusionForce<2> diffusion_force;
        TS_ASSERT_EQUALS(diffusion_force.GetIdentifier(), "DiffusionForce-2");

        out_stream diffusion_force_parameter_file = output_file_handler.OpenOutputFile("diffusion_results.parameters");
        diffusion_force.OutputForceParameters(diffusion_force_parameter_file);
        diffusion_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("diffusion_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/diffusion_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestGeneralisedLinearSpringForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
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

    void TestDifferentialAdhesionGeneralisedLinearSpringForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DifferentialAdhesionGeneralisedLinearSpringForce.arch";

        {
            DifferentialAdhesionGeneralisedLinearSpringForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetMeinekeSpringStiffness(12.34);
            force.SetMeinekeDivisionRestingSpringLength(0.856);
            force.SetMeinekeSpringGrowthDuration(2.593);
            force.SetHomotypicLabelledSpringConstantMultiplier(0.051);
            force.SetHeterotypicSpringConstantMultiplier(1.348);

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
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringStiffness(), 12.34, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionGeneralisedLinearSpringForce<2>*>(p_force))->GetHomotypicLabelledSpringConstantMultiplier(), 0.051, 1e-6);
            TS_ASSERT_DELTA((static_cast<DifferentialAdhesionGeneralisedLinearSpringForce<2>*>(p_force))->GetHeterotypicSpringConstantMultiplier(), 1.348, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestChemotacticForceMethods()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        unsigned cells_across = 7;
        unsigned cells_up = 5;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        HoneycombMeshGenerator generator(cells_across, cells_up);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MAKE_PTR(CellLabel, p_label);
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-10);
            cells[i]->AddCellProperty(p_label);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell data on the cell population
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(p_mesh->GetNode(i)->GetIndex());
            p_cell->GetCellData()->SetItem("nutrient", x/50.0);

        }

        ChemotacticForce<2> chemotactic_force;

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        chemotactic_force.AddForceContribution(cell_population);

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
            TS_ASSERT_DELTA(cell_population.GetNode(index)->rGetAppliedForce()[0], force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(index)->rGetAppliedForce()[1], 0.0, 1e-4);
        }
    }

    void TestChemotacticForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
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

    void TestRepulsionForceMethods()
    {
        // Create a NodeBasedCellPopulation
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(2, true, 3.0, 0.0));

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 100.0);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Update(); //Needs to be called separately as not in a simulation

        RepulsionForce<2> repulsion_force;

        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }
        repulsion_force.AddForceContribution(cell_population);

        /*
         * First two cells repel each other and second 2 cells are too far apart.
         * The radius of the cells is the default value, 0.5.
         */
        if (PetscTools::AmMaster())    // All cells in this test lie on the master process.
        {
            unsigned zero_index = 0;
            unsigned one_index = PetscTools::GetNumProcs();
            unsigned two_index = 2*PetscTools::GetNumProcs();
            TS_ASSERT_DELTA(cell_population.GetNode(zero_index)->rGetAppliedForce()[0], -34.5387, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(zero_index)->rGetAppliedForce()[1], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(one_index)->rGetAppliedForce()[0], 34.5387, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(one_index)->rGetAppliedForce()[1], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(two_index)->rGetAppliedForce()[0], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(two_index)->rGetAppliedForce()[1], 0.0, 1e-4);

            // Tests the calculation of the force with different cell radii
            mesh.GetNode(zero_index)->SetRadius(10);
            mesh.GetNode(one_index)->SetRadius(10);
            mesh.GetNode(two_index)->SetRadius(10);

            // Reset the vector of node forces
            for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                    node_iter != mesh.GetNodeIteratorEnd();
                    ++node_iter)
            {
                node_iter->ClearAppliedForce();
            }

            repulsion_force.AddForceContribution(cell_population);

            // All cells repel each other
            TS_ASSERT_DELTA(cell_population.GetNode(zero_index)->rGetAppliedForce()[0], 15.0 * 20.0 * log(1 - 19.9/20)+15.0 * 20.0 * log(1 - 17.0/20), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(zero_index)->rGetAppliedForce()[1], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(one_index)->rGetAppliedForce()[0], -15.0 * 20.0 * log(1 - 19.9/20)+15.0 * 20.0 * log(1 - 17.1/20), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(one_index)->rGetAppliedForce()[1], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(two_index)->rGetAppliedForce()[0], -15.0 * 20.0 * log(1 - 17.1/20)-15.0 * 20.0 * log(1 - 17.0/20), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(two_index)->rGetAppliedForce()[1], 0.0, 1e-4);

            // Tests the calculation of the force with different cell radii
            mesh.GetNode(zero_index)->SetRadius(0.2);
            mesh.GetNode(one_index)->SetRadius(0.2);
            mesh.GetNode(two_index)->SetRadius(0.2);

            // Reset the vector of node forces
            for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                    node_iter != mesh.GetNodeIteratorEnd();
                    ++node_iter)
            {
                node_iter->ClearAppliedForce();
            }

            repulsion_force.AddForceContribution(cell_population);

            /*
             * First two cells repel each other and second 2 cells are too far apart.
             * The overlap is -0.3 and the spring stiffness is the default value 15.0.
             */
            TS_ASSERT_DELTA(cell_population.GetNode(zero_index)->rGetAppliedForce()[0], 15.0 * 0.4 * log(1 -0.3/0.4), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(zero_index)->rGetAppliedForce()[1], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(one_index)->rGetAppliedForce()[0], -15.0 * 0.4 * log(1 -0.3/0.4), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(one_index)->rGetAppliedForce()[1], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(two_index)->rGetAppliedForce()[0], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(two_index)->rGetAppliedForce()[1], 0.0, 1e-4);
        }

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestRepulsionForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
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

    void TestNagaiHondaForceMethods()
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
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
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
        TS_ASSERT_DELTA(force.GetNagaiHondaCellCellAdhesionEnergyParameter(), 0.5, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaCellBoundaryAdhesionEnergyParameter(), 1.0, 1e-12);

        force.SetNagaiHondaDeformationEnergyParameter(5.8);
        force.SetNagaiHondaMembraneSurfaceEnergyParameter(17.9);
        force.SetNagaiHondaCellCellAdhesionEnergyParameter(0.3);
        force.SetNagaiHondaCellBoundaryAdhesionEnergyParameter(0.6);

        TS_ASSERT_DELTA(force.GetNagaiHondaDeformationEnergyParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaMembraneSurfaceEnergyParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaCellCellAdhesionEnergyParameter(), 0.3, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-12);

        force.SetNagaiHondaDeformationEnergyParameter(100.0);
        force.SetNagaiHondaMembraneSurfaceEnergyParameter(10.0);
        force.SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        force.SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Currently, NagaiHonda force only works if used together with a target area growth modifier
        // This tests that a meaningful error appears if we don't use a growth modifier
        TS_ASSERT_THROWS_THIS(force.AddForceContribution(cell_population),
                "You need to add an AbstractTargetAreaModifier to the simulation in order to use NagaiHondaForce");

        // create our modifier, which sets the target areas for the cell population
        // this is a workaround to make the test work
        // #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should be radially inward, with the same magnitude for all nodes
        double force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.25, 2);

        // Set the cell to be necrotic
        cell_population.GetCellUsingLocationIndex(0)->StartApoptosis();

        // Reset force vector
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should not yet be affected by setting the cell to be apoptotic
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Increment time
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_DELTA(cell_population.GetCellUsingLocationIndex(0)->GetTimeUntilDeath(), 0.125, 1e-6);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // Now the forces should be affected
        double apoptotic_force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());
        TS_ASSERT_LESS_THAN(force_magnitude, apoptotic_force_magnitude);
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), apoptotic_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -apoptotic_force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -apoptotic_force_magnitude*sin(angles[i]), 1e-4);
        }
    }

    void TestNagaiHondaForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
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

            // Tidy up
            delete p_force;
        }
    }

    void TestNagaiHondaDifferentialAdhesionForceMethods()
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        MAKE_PTR(CellLabel, p_label);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force system
        NagaiHondaDifferentialAdhesionForce<2> force;

        // Set member variables
        force.SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(9.1);
        force.SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(2.8);
        force.SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(7.3);
        force.SetNagaiHondaCellCellAdhesionEnergyParameter(6.4);
        force.SetNagaiHondaCellBoundaryAdhesionEnergyParameter(0.6);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // now we add the the growth modifier and go on
        // #2488
        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // Test the case where the nodes are shared by a cell on the boundary that is not labelled
        Node<2>* p_node_0 = p_mesh->GetElement(0)->GetNode(0);
        Node<2>* p_node_4 = p_mesh->GetElement(0)->GetNode(1);
        double adhesion_parameter_nodes_0_4 = force.GetAdhesionParameter(p_node_0, p_node_4, cell_population);
        TS_ASSERT_DELTA(adhesion_parameter_nodes_0_4, 0.6, 1e-6);

        // Test the case where the nodes are shared by 3 non-labelled cells
        Node<2>* p_node_10 = p_mesh->GetElement(4)->GetNode(0);
        Node<2>* p_node_14 = p_mesh->GetElement(4)->GetNode(1);
        double adhesion_parameter_nodes_10_14 = force.GetAdhesionParameter(p_node_10, p_node_14, cell_population);
        TS_ASSERT_DELTA(adhesion_parameter_nodes_10_14, 6.4, 1e-6);

        // Test the case where the nodes are shared by a cell on the boundary that is labelled
        cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        adhesion_parameter_nodes_0_4 = force.GetAdhesionParameter(p_node_0, p_node_4, cell_population);
        TS_ASSERT_DELTA(adhesion_parameter_nodes_0_4, 7.3, 1e-6);

        // Test the case where the nodes are shared by 1 labelled cell
        cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_label);
        adhesion_parameter_nodes_10_14 = force.GetAdhesionParameter(p_node_10, p_node_14, cell_population);
        TS_ASSERT_DELTA(adhesion_parameter_nodes_10_14, 9.1, 1e-6);

        // Test the case where the nodes are shared by 2 labelled cells
        for (unsigned i=0; i<cell_population.GetNumElements(); i++)
        {
            cell_population.GetCellUsingLocationIndex(i)->AddCellProperty(p_label);
        }
        adhesion_parameter_nodes_10_14 = force.GetAdhesionParameter(p_node_10, p_node_14, cell_population);
        TS_ASSERT_DELTA(adhesion_parameter_nodes_10_14, 2.8, 1e-6);
    }

    void TestNagaiHondaDifferentialAdhesionForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "NagaiHondaDifferentialAdhesionForce.arch";

        {
            NagaiHondaDifferentialAdhesionForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(9.1);
            force.SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(2.8);
            force.SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(7.3);

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
            TS_ASSERT_DELTA(static_cast<NagaiHondaDifferentialAdhesionForce<2>*>(p_force)->GetNagaiHondaLabelledCellCellAdhesionEnergyParameter(), 9.1, 1e-12);
            TS_ASSERT_DELTA(static_cast<NagaiHondaDifferentialAdhesionForce<2>*>(p_force)->GetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(), 2.8, 1e-12);
            TS_ASSERT_DELTA(static_cast<NagaiHondaDifferentialAdhesionForce<2>*>(p_force)->GetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(), 7.3, 1e-12);

            // Tidy up
            delete p_force;
        }
    }

    void TestWelikyOsterForceMethods()
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
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
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

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        force.AddForceContribution(cell_population);

        // The force on each node should be radially inward, with the same magnitude for all nodes
        double force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }
    }

    void TestWelikyOsterForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
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

    void TestFarhadifarForceMethods()
    {
        // This is the same test as for other vertex based forces. It comprises a sanity check that forces point in the right direction.
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
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(-1.0);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.InitialiseCells();

        // Create a force system
        FarhadifarForce<2> force;

        // Test get/set methods
        TS_ASSERT_DELTA(force.GetAreaElasticityParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(force.GetPerimeterContractilityParameter(), 0.04, 1e-12);
        TS_ASSERT_DELTA(force.GetLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetBoundaryLineTensionParameter(), 0.12, 1e-12);

        force.SetAreaElasticityParameter(5.8);
        force.SetPerimeterContractilityParameter(17.9);
        force.SetLineTensionParameter(0.5);
        force.SetBoundaryLineTensionParameter(0.6);

        TS_ASSERT_DELTA(force.GetAreaElasticityParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(force.GetPerimeterContractilityParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(force.GetLineTensionParameter(), 0.5, 1e-12);
        TS_ASSERT_DELTA(force.GetBoundaryLineTensionParameter(), 0.6, 1e-12);

        force.SetAreaElasticityParameter(1.0);
        force.SetPerimeterContractilityParameter(0.04);
        force.SetLineTensionParameter(0.12);
        force.SetBoundaryLineTensionParameter(0.12);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Currently, the Farhadifar force only works if used together with a target area growth modifier
        // This tests that a meaningful error appears if we don't use a growth modifier
        TS_ASSERT_THROWS_THIS(force.AddForceContribution(cell_population),
                "You need to add an AbstractTargetAreaModifier to the simulation in order to use a FarhadifarForce");

        // create our modifier, which sets the target areas for the cell population

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should be radially inward, with the same magnitude for all nodes
        double force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.25, 2);

        // Set the cell to be necrotic
        cell_population.GetCellUsingLocationIndex(0)->StartApoptosis();

        // Reset force vector
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should not yet be affected by setting the cell to be apoptotic
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Increment time
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_DELTA(cell_population.GetCellUsingLocationIndex(0)->GetTimeUntilDeath(), 0.125, 1e-6);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // Now the forces should be affected
        double apoptotic_force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());
        TS_ASSERT_LESS_THAN(force_magnitude, apoptotic_force_magnitude);
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), apoptotic_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -apoptotic_force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -apoptotic_force_magnitude*sin(angles[i]), 1e-4);
        }
    }

    void TestFarhadifarForceTerms()
       {
        /**
         * Here we test that the forces are applied correctly to individual nodes.
         * We apply the force to something like this:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
        }

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        // Now let's make a FarhadifarForce and apply it to the population
        FarhadifarForce<2> force;

        force.AddForceContribution(cell_population);

        c_vector<double, 2> applied_force_0;
        applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1;
        applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // If this is a Farhadifar force, this will be the force at the vertices
        TS_ASSERT_DELTA(applied_force_0[0], 3.44, 1e-10);
        TS_ASSERT_DELTA(applied_force_0[1], 3.44, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[0], 0.0, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[1], 6.76, 1e-10);
    }

    void TestFarhadifarForceInSimulation()
    {
        /**
         * This is the same test as above, just that now we don't check that the applied forces are calculated correctly,
         * but rather that in a simulation the displacement of vertices is as we expect.
         *
         * This is the mesh:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
        }

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        // Now let's make a FarhadifarForce and add it to the simulation.
        MAKE_PTR(FarhadifarForce<2>, p_force);

        // We need to reset the cell rearrangement threshold - vertex movements are kept below that threshold
        cell_population.rGetMesh().SetCellRearrangementThreshold(0.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFarhadifarForce");
        simulator.SetEndTime(0.01);
        simulator.SetDt(0.01);
        simulator.AddForce(p_force);

        simulator.Solve();

        c_vector<double, 2> applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // New Location = Old Location + (Dt * applied force), since viscosity should be one
        c_vector<double, 2> expected_new_node_location_0;
        expected_new_node_location_0[0] = 0.0+0.01*3.44;
        expected_new_node_location_0[1] = 0.0+0.01*3.44;
        c_vector<double, 2> expected_new_node_location_1;
        expected_new_node_location_1[0] = 2.0 + 0.01*0.0;
        expected_new_node_location_1[1] = 0.0 + 0.01*6.76;

        // If this is a Farhadifar force, this will be the location of the first two vertices.
        TS_ASSERT_DELTA(expected_new_node_location_0[0], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_0[1], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[1], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[0], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[1], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[1], 1e-10);

    }

    void TestFarhadifarForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "FarhadifarForce.arch";

        {
            FarhadifarForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetAreaElasticityParameter(5.8);
            force.SetPerimeterContractilityParameter(17.9);
            force.SetLineTensionParameter(0.5);
            force.SetBoundaryLineTensionParameter(0.6);

            // Serialize via pointer to most abstract class possible
            AbstractForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<2>* p_abstract_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_abstract_force;

            FarhadifarForce<2>* p_farhadifar_force = static_cast<FarhadifarForce<2>*>(p_abstract_force);

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(p_farhadifar_force->GetAreaElasticityParameter(), 5.8, 1e-12);
            TS_ASSERT_DELTA(p_farhadifar_force->GetPerimeterContractilityParameter(), 17.9, 1e-12);
            TS_ASSERT_DELTA(p_farhadifar_force->GetLineTensionParameter(), 0.5, 1e-12);
            TS_ASSERT_DELTA(p_farhadifar_force->GetBoundaryLineTensionParameter(), 0.6, 1e-12);

            // Tidy up
            delete p_abstract_force;
        }
    }

    void TestCentreBasedForcesWithVertexCellPopulation()
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
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(-1.0);
        cells.push_back(p_cell);

        // Create VertexBasedCellPopulation
        VertexBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.InitialiseCells();

        // Test that a subclass of AbstractTwoBodyInteractionForce throws the correct exception
        GeneralisedLinearSpringForce<2> spring_force;
        TS_ASSERT_THROWS_THIS(spring_force.AddForceContribution(cell_population),
                "Subclasses of AbstractTwoBodyInteractionForce are to be used with subclasses of AbstractCentreBasedCellPopulation only");

        // Test that RepulsionForce throws the correct exception
        RepulsionForce<2> repulsion_force;
        TS_ASSERT_THROWS_THIS(repulsion_force.AddForceContribution(cell_population),
                 "RepulsionForce is to be used with a NodeBasedCellPopulation only");
    }

    void TestIncorrectForcesWithNodeBasedCellPopulation()
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
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Test that NagaiHondaForce throws the correct exception
        NagaiHondaForce<2> nagai_honda_force;
        TS_ASSERT_THROWS_THIS(nagai_honda_force.AddForceContribution(cell_population),
                "NagaiHondaForce is to be used with a VertexBasedCellPopulation only");

        // Test that WelikyOsterForce throws the correct exception
        WelikyOsterForce<2> weliky_oster_force;
        TS_ASSERT_THROWS_THIS(weliky_oster_force.AddForceContribution(cell_population),
                "WelikyOsterForce is to be used with a VertexBasedCellPopulation only");

        // Test that FarhadifarForce throws the correct exception
        FarhadifarForce<2> farhadifar_force;
        TS_ASSERT_THROWS_THIS(farhadifar_force.AddForceContribution(cell_population),
                "FarhadifarForce is to be used with a VertexBasedCellPopulation only");

        // Avoid memory leak
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
        MutableMesh<1,1> generating_mesh;
        generating_mesh.ConstructLinearMesh(5);

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes(), std::vector<unsigned>(), p_diff_type);

        // Create cell population
        std::vector<CellPtr> cells_copy(cells);
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Create force law object
        DiffusionForce<1> diffusion_force;

        for (AbstractMesh<1,1>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        // Compute forces on nodes
        diffusion_force.AddForceContribution(cell_population);

        // Test Set and Get methods for the diffusion force
        TS_ASSERT_DELTA(diffusion_force.GetViscosity(), 3.204e-6, 1e-10);
        TS_ASSERT_DELTA(diffusion_force.GetAbsoluteTemperature(), 296.0, 1e-10);

        diffusion_force.SetViscosity(0.01);
        diffusion_force.SetAbsoluteTemperature(100.0);
        TS_ASSERT_DELTA(diffusion_force.GetViscosity(), 0.01, 1e-10);
        TS_ASSERT_DELTA(diffusion_force.GetAbsoluteTemperature(), 100.0, 1e-10);
        diffusion_force.SetViscosity(3.204e-6);
        diffusion_force.SetAbsoluteTemperature(296.0);
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
        mesh.ConstructNodesWithoutMesh(nodes, 100.0);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Update(); //Needs to be called separately as not in a simulation

        // Create force law object
        DiffusionForce<2> force;

        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        if (mesh.GetNumNodes() > 0)
        {
            // Loop over time iterations
            double variance = 0.0;
            unsigned num_iterations = 1000;
            for (unsigned i=0; i<num_iterations; i++)
            {
                // Re-initialize the force on node zero
                mesh.GetNodeIteratorBegin()->ClearAppliedForce();

                // Compute forces on nodes
                force.AddForceContribution(cell_population);

                // Calculate the variance
                variance += pow(norm_2(mesh.GetNodeIteratorBegin()->rGetAppliedForce()),2);
            }

            double correct_diffusion_coefficient =
                    4.97033568e-7 * force.GetAbsoluteTemperature() / (6 * M_PI * force.GetViscosity() * mesh.GetNodeIteratorBegin()->GetRadius() );
            unsigned dim = 2;
            variance /= num_iterations*2*dim*correct_diffusion_coefficient*SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(variance, 1.0, 1e-1);
        }

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

        // Tidy up
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    void TestDiffusionForceWithVertexBasedCellPopulation()
    {
        // Define the seed
        RandomNumberGenerator::Instance()->Reseed(0);

        // Set up time parameters
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create a simple VertexBasedCellPopulation
        HoneycombVertexMeshGenerator mesh_generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();
        for (AbstractMesh<2,2>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
             node_iter != p_mesh->GetNodeIteratorEnd();
             ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create DiffusionForce object
        DiffusionForce<2> force;

        // Check that AddForceContribution() throws the right error if the node radii have not been set
        TS_ASSERT_THROWS_THIS(force.AddForceContribution(cell_population),
            "SetRadius() must be called on each Node before calling DiffusionForce::AddForceContribution() to avoid a division by zero error");

        // Now set each node radius...
        for (AbstractMesh<2,2>::NodeIterator node_iter = cell_population.rGetMesh().GetNodeIteratorBegin();
             node_iter != cell_population.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            node_iter->SetRadius(1.0);
        }

        // ...and check that AddForceContribution() throws no error
        TS_ASSERT_THROWS_NOTHING(force.AddForceContribution(cell_population));

        // Tidy up
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    void TestDiffusionForceWithMeshBasedCellPopulation()
    {
        EXIT_IF_PARALLEL;

        // Define the seed
        RandomNumberGenerator::Instance()->Reseed(0);

        // Set up time parameters
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create a simple MeshBasedCellPopulation
        HoneycombMeshGenerator mesh_generator(4, 6, 0);
        MutableMesh<2,2>* p_mesh = mesh_generator.GetMesh();
        for (AbstractMesh<2,2>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
             node_iter != p_mesh->GetNodeIteratorEnd();
             ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes(), std::vector<unsigned>(), p_diff_type);
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create DiffusionForce object
        DiffusionForce<2> force;

        // Check that AddForceContribution() throws the right error if the node radii have not been set
        TS_ASSERT_THROWS_THIS(force.AddForceContribution(cell_population),
            "SetRadius() must be called on each Node before calling DiffusionForce::AddForceContribution() to avoid a division by zero error");

        // Now set each node radius...
        for (AbstractMesh<2,2>::NodeIterator node_iter = cell_population.rGetMesh().GetNodeIteratorBegin();
             node_iter != cell_population.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            node_iter->SetRadius(1.0);
        }

        // ...and check that AddForceContribution() throws no error
        TS_ASSERT_THROWS_NOTHING(force.AddForceContribution(cell_population));

        // Now update the cell population, which in turn calls ReMesh() on the mesh...
        cell_population.Update();

        for (AbstractMesh<2,2>::NodeIterator node_iter = cell_population.rGetMesh().GetNodeIteratorBegin();
             node_iter != cell_population.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            TS_ASSERT_DELTA(node_iter->GetRadius(), 1.0, 1e-6);
        }

        // ...and check that AddForceContribution() still throws no error
        TS_ASSERT_THROWS_NOTHING(force.AddForceContribution(cell_population));

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
        mesh.ConstructNodesWithoutMesh(nodes, 100.0);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.Update(); //Needs to be called separately as not in a simulation

        // Create force law object
        DiffusionForce<3> force;

        for (AbstractMesh<3,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        double variance = 0.0;

        // Loop over time iterations
        if (mesh.GetNumNodes() > 0)
        {
            unsigned num_iterations = 1000;
            for (unsigned i=0; i<num_iterations; i++)
            {
                // Re-initialize the force on node zero
                mesh.GetNodeIteratorBegin()->ClearAppliedForce();

                // Compute forces on nodes
                force.AddForceContribution(cell_population);

                // Calculate the variance
                variance += pow(norm_2(mesh.GetNodeIteratorBegin()->rGetAppliedForce()),2);
            }

            double correct_diffusion_coefficient =
                    4.97033568e-7 * force.GetAbsoluteTemperature() / (6 * M_PI * force.GetViscosity() * mesh.GetNodeIteratorBegin()->GetRadius() );
            unsigned dim = 3;
            variance /= num_iterations*2*dim*correct_diffusion_coefficient*SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(variance, 1.0, 1e-1);
        }

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

        // Tidy up
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    void TestDiffusionForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
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
            TS_ASSERT_DELTA((static_cast<DiffusionForce<2>*>(p_force))->GetAbsoluteTemperature(), 296.0, 1e-6);
            TS_ASSERT_DELTA((static_cast<DiffusionForce<2>*>(p_force))->GetViscosity(), 3.204e-6, 1e-6);

            // Tidy up
            delete p_force;
        }
    }
};

#endif /*TESTFORCES_HPP_*/
