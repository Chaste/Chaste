/*

Copyright (c) 2005-2026, University of Oxford.
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
#include "PlanarPolarisedFarhadifarForce.hpp"
#include "DiffusionForce.hpp"
#include "SemForce.hpp"
#include "SemLinearForce.hpp"
#include "SemGaussianRandomForce.hpp"
#include "SemSpatiallyCorrelatedRandomForce.hpp"
#include "SemSingleElementMeshGenerator.hpp"
#include "SemMultiElementMeshGenerator.hpp"
#include "SemBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
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
#include "Warnings.hpp"

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
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();
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
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();
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

        // Test with FarhadifarForce
        FarhadifarForce<2> farhadifar_force;
        TS_ASSERT_EQUALS(farhadifar_force.GetIdentifier(), "FarhadifarForce-2");

        out_stream farhadifar_force_parameter_file = output_file_handler.OpenOutputFile("farhadifar_results.parameters");
        farhadifar_force.OutputForceParameters(farhadifar_force_parameter_file);
        farhadifar_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("farhadifar_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/farhadifar_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with PlanarPolarisedFarhadifarForce
        PlanarPolarisedFarhadifarForce<2> planar_force;
        TS_ASSERT_EQUALS(planar_force.GetIdentifier(), "PlanarPolarisedFarhadifarForce-2");

        out_stream planar_force_parameter_file = output_file_handler.OpenOutputFile("planar_results.parameters");
        planar_force.OutputForceParameters(planar_force_parameter_file);
        planar_force_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("planar_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestForces/planar_results.parameters",
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
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

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
        TS_ASSERT_DELTA(force.GetNagaiHondaTargetAreaParameter(), 1.0, 1e-12);

        force.SetNagaiHondaDeformationEnergyParameter(5.8);
        force.SetNagaiHondaMembraneSurfaceEnergyParameter(17.9);
        force.SetNagaiHondaCellCellAdhesionEnergyParameter(0.3);
        force.SetNagaiHondaCellBoundaryAdhesionEnergyParameter(0.6);
        force.SetNagaiHondaTargetAreaParameter(1.7);

        TS_ASSERT_DELTA(force.GetNagaiHondaDeformationEnergyParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaMembraneSurfaceEnergyParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaCellCellAdhesionEnergyParameter(), 0.3, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaCellBoundaryAdhesionEnergyParameter(), 0.6, 1e-12);
        TS_ASSERT_DELTA(force.GetNagaiHondaTargetAreaParameter(), 1.7, 1e-12);

        force.SetNagaiHondaDeformationEnergyParameter(100.0);
        force.SetNagaiHondaMembraneSurfaceEnergyParameter(10.0);
        force.SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        force.SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
        force.SetNagaiHondaTargetAreaParameter(1.0);

        // Calculate force on each node
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
        force.AddForceContribution(cell_population);

        // The force on each node should not yet be affected by setting the cell to be apoptotic
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Modify cell target areas over time according to a simple growth model
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

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
            force.SetNagaiHondaTargetAreaParameter(3.2);

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
            TS_ASSERT_DELTA(static_cast<NagaiHondaForce<2>*>(p_force)->GetNagaiHondaTargetAreaParameter(), 3.2, 1e-12);

            // Tidy up
            delete p_force;
        }
    }

    void TestNagaiHondaDifferentialAdhesionForceMethods()
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(3, 3);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

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
        TS_ASSERT_DELTA(force.GetTargetAreaParameter(), 1.0, 1e-12);

        force.SetAreaElasticityParameter(5.8);
        force.SetPerimeterContractilityParameter(17.9);
        force.SetLineTensionParameter(0.5);
        force.SetBoundaryLineTensionParameter(0.6);
        force.SetTargetAreaParameter(2.9);

        TS_ASSERT_DELTA(force.GetAreaElasticityParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(force.GetPerimeterContractilityParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(force.GetLineTensionParameter(), 0.5, 1e-12);
        TS_ASSERT_DELTA(force.GetBoundaryLineTensionParameter(), 0.6, 1e-12);
        TS_ASSERT_DELTA(force.GetTargetAreaParameter(), 2.9, 1e-12);

        force.SetAreaElasticityParameter(1.0);
        force.SetPerimeterContractilityParameter(0.04);
        force.SetLineTensionParameter(0.12);
        force.SetBoundaryLineTensionParameter(0.12);
        force.SetTargetAreaParameter(1.0);

        // Calculate force on each node
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

        force.AddForceContribution(cell_population);

        // The force on each node should not yet be affected by setting the cell to be apoptotic
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Modify cell target areas over time according to a simple growth model
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

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

    void TestPlanarPolarisedFarhadifarForce()
    {
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 9;
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

        // Create force
        PlanarPolarisedFarhadifarForce<2> force;

        // Test get/set methods
        TS_ASSERT_DELTA(force.GetAreaElasticityParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(force.GetPerimeterContractilityParameter(), 0.04, 1e-12);
        TS_ASSERT_DELTA(force.GetLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetBoundaryLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetPlanarPolarisedLineTensionMultiplier(), 2.0, 1e-12);

        force.SetPlanarPolarisedLineTensionMultiplier(4.0);
        TS_ASSERT_DELTA(force.GetPlanarPolarisedLineTensionMultiplier(), 4.0, 1e-12);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Test GetLineTensionParameter()
        c_vector<double, 2> applied_force_0;
        applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1;
        applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        for (unsigned node_idx = 0; node_idx < cell_population.GetNumNodes(); node_idx++)
        {
            Node<2>* p_node_A = cell_population.GetNode(node_idx);
            Node<2>* p_node_B = cell_population.GetNode((node_idx + 1) % cell_population.GetNumNodes());

            double line_tension = force.GetLineTensionParameter(p_node_A, p_node_B, cell_population);
            if ((node_idx == 1) || (node_idx == 2) || (node_idx == 6) || (node_idx == 7))
            {
                TS_ASSERT_DELTA(0.12, line_tension, 1e-12);
            }
            else
            {
                TS_ASSERT_DELTA(0.48, line_tension, 1e-12);
            }
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

    void TestPlanarPolarisedFarhadifarForceArchiving()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "PlanarPolarisedFarhadifarForce.arch";

        {
            PlanarPolarisedFarhadifarForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetAreaElasticityParameter(5.8);
            force.SetPerimeterContractilityParameter(17.9);
            force.SetLineTensionParameter(0.5);
            force.SetBoundaryLineTensionParameter(0.6);
            force.SetPlanarPolarisedLineTensionMultiplier(5.2);

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

            PlanarPolarisedFarhadifarForce<2>* p_planar_force = static_cast<PlanarPolarisedFarhadifarForce<2>*>(p_abstract_force);

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(p_planar_force->GetAreaElasticityParameter(), 5.8, 1e-12);
            TS_ASSERT_DELTA(p_planar_force->GetPerimeterContractilityParameter(), 17.9, 1e-12);
            TS_ASSERT_DELTA(p_planar_force->GetLineTensionParameter(), 0.5, 1e-12);
            TS_ASSERT_DELTA(p_planar_force->GetBoundaryLineTensionParameter(), 0.6, 1e-12);
            TS_ASSERT_DELTA(p_planar_force->GetPlanarPolarisedLineTensionMultiplier(), 5.2, 1e-12);

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

        // Test that a subclass of AbstractTwoBodyInteractionForce emits the appropriate warning
        GeneralisedLinearSpringForce<2> spring_force;
        Warnings::QuietDestroy(); // Clear any warnings before the expected one
        spring_force.AddForceContribution(cell_population);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),
            "No node pairs found. Does this cell population support updating mNodePairs?"
            " Currently this force class works only with NodeBased and MeshBased cell populations.");

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
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = mesh_generator.GetMesh();
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
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = mesh_generator.GetMesh();
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

    /**
     * Test getter and setter methods for SemForce.
     */
    void TestSemForceGettersSetters()
    {
        SemForce<2> force;

        // Test default values
        TS_ASSERT_DELTA(force.GetIntraWellDepth(), 1.0, 1e-6);
        TS_ASSERT_DELTA(force.GetIntraScalingFactor(), 5.0, 1e-6);
        TS_ASSERT_DELTA(force.GetIntraEquilibriumDistance(), 0.2, 1e-6);
        TS_ASSERT_DELTA(force.GetIntraCutOffDistance(), 0.5, 1e-6);
        TS_ASSERT_DELTA(force.GetInterWellDepth(), 1.0, 1e-6);
        TS_ASSERT_DELTA(force.GetInterScalingFactor(), 5.0, 1e-6);
        TS_ASSERT_DELTA(force.GetInterEquilibriumDistance(), 0.3, 1e-6);
        TS_ASSERT_DELTA(force.GetInterCutOffDistance(), 0.5, 1e-6);

        // Test setters
        force.SetIntraWellDepth(2.0);
        force.SetIntraScalingFactor(3.0);
        force.SetIntraEquilibriumDistance(0.15);
        force.SetIntraCutOffDistance(0.6);
        force.SetInterWellDepth(1.5);
        force.SetInterScalingFactor(4.0);
        force.SetInterEquilibriumDistance(0.25);
        force.SetInterCutOffDistance(0.7);

        TS_ASSERT_DELTA(force.GetIntraWellDepth(), 2.0, 1e-6);
        TS_ASSERT_DELTA(force.GetIntraScalingFactor(), 3.0, 1e-6);
        TS_ASSERT_DELTA(force.GetIntraEquilibriumDistance(), 0.15, 1e-6);
        TS_ASSERT_DELTA(force.GetIntraCutOffDistance(), 0.6, 1e-6);
        TS_ASSERT_DELTA(force.GetInterWellDepth(), 1.5, 1e-6);
        TS_ASSERT_DELTA(force.GetInterScalingFactor(), 4.0, 1e-6);
        TS_ASSERT_DELTA(force.GetInterEquilibriumDistance(), 0.25, 1e-6);
        TS_ASSERT_DELTA(force.GetInterCutOffDistance(), 0.7, 1e-6);
    }

    /**
     * Test the Morse force calculation directly via CalculateForceVector.
     * Verifies equilibrium, repulsive/attractive sign convention, and exact values.
     */
    void TestSemMorseForceCalculation()
    {
        SemForce<2> morse_force;

        double u0 = 1.0;
        double rho = 5.0;
        double r_eq = 0.2;

        // Test 1: At equilibrium distance, force should be zero
        {
            c_vector<double, 2> vec_a_to_b;
            vec_a_to_b[0] = r_eq;
            vec_a_to_b[1] = 0.0;
            double dist_sq = r_eq * r_eq;

            c_vector<double, 2> force = morse_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            TS_ASSERT_DELTA(force[0], 0.0, 1e-10);
            TS_ASSERT_DELTA(force[1], 0.0, 1e-10);
        }

        // Test 2: At r < r_eq, force should be repulsive (pointing away from B, i.e. negative x)
        {
            double r = 0.15;
            c_vector<double, 2> vec_a_to_b;
            vec_a_to_b[0] = r;
            vec_a_to_b[1] = 0.0;
            double dist_sq = r * r;

            c_vector<double, 2> force = morse_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            TS_ASSERT_LESS_THAN(force[0], 0.0); // Repulsive: force on A away from B
            TS_ASSERT_DELTA(force[1], 0.0, 1e-10);
        }

        // Test 3: At r > r_eq, force should be attractive (pointing toward B, i.e. positive x)
        {
            double r = 0.25;
            c_vector<double, 2> vec_a_to_b;
            vec_a_to_b[0] = r;
            vec_a_to_b[1] = 0.0;
            double dist_sq = r * r;

            c_vector<double, 2> force = morse_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            TS_ASSERT_LESS_THAN(0.0, force[0]); // Attractive: force on A toward B
            TS_ASSERT_DELTA(force[1], 0.0, 1e-10);
        }

        // Test 4: Verify exact value against hand calculation
        // At r = 0.15: s = 0.15^2/0.2^2 = 0.5625
        // exp(rho*(1-s)) = exp(5*0.4375) = exp(2.1875)
        // coeff = (4*5*1/0.04) * (exp(2.1875) - exp(2*2.1875))
        //       = 500 * (exp(2.1875) - exp(4.375))
        // force_x = coeff * 0.15
        {
            double r = 0.15;
            c_vector<double, 2> vec_a_to_b;
            vec_a_to_b[0] = r;
            vec_a_to_b[1] = 0.0;
            double dist_sq = r * r;

            double s = dist_sq / (r_eq * r_eq);
            double exp_rho_val = std::exp(rho * (1.0 - s));
            double exp_2rho_val = exp_rho_val * exp_rho_val;
            double expected_coeff = (4.0 * rho * u0 / (r_eq * r_eq)) * (exp_rho_val - exp_2rho_val);
            double expected_force_x = expected_coeff * r;

            c_vector<double, 2> force = morse_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            TS_ASSERT_DELTA(force[0], expected_force_x, 1e-6);
            TS_ASSERT_DELTA(force[1], 0.0, 1e-10);
        }

        // Test 5: Force should be along the vector A->B in a diagonal direction
        {
            double r = 0.25;
            c_vector<double, 2> vec_a_to_b;
            vec_a_to_b[0] = r / std::sqrt(2.0);
            vec_a_to_b[1] = r / std::sqrt(2.0);
            double dist_sq = r * r;

            c_vector<double, 2> force = morse_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            // Force should be symmetric in x and y
            TS_ASSERT_DELTA(force[0], force[1], 1e-10);
            // And attractive (positive, toward B)
            TS_ASSERT_LESS_THAN(0.0, force[0]);
        }
    }

    /**
     * Test the linear (harmonic) force calculation and verify it agrees
     * with the full Morse force for small deviations about equilibrium.
     */
    void TestSemLinearForceCalculation()
    {
        SemLinearForce<2> linear_force;

        double u0 = 1.0;
        double rho = 5.0;
        double r_eq = 0.2;
        double kappa = 8.0 * rho * rho * u0 / (r_eq * r_eq);

        // Test 1: At equilibrium, force should be zero
        {
            c_vector<double, 2> vec_a_to_b;
            vec_a_to_b[0] = r_eq;
            vec_a_to_b[1] = 0.0;
            double dist_sq = r_eq * r_eq;

            c_vector<double, 2> force = linear_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            TS_ASSERT_DELTA(force[0], 0.0, 1e-10);
            TS_ASSERT_DELTA(force[1], 0.0, 1e-10);
        }

        // Test 2: Verify exact value at r = 0.25 (stretched by 0.05)
        // F_A = kappa * (1 - r_eq/r) * vec_a_to_b
        //     = kappa * (1 - 0.2/0.25) * 0.25 = kappa * 0.2 * 0.25 = kappa * 0.05
        {
            double r = 0.25;
            c_vector<double, 2> vec_a_to_b;
            vec_a_to_b[0] = r;
            vec_a_to_b[1] = 0.0;
            double dist_sq = r * r;

            double expected_force_x = kappa * (1.0 - r_eq / r) * r;

            c_vector<double, 2> force = linear_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            TS_ASSERT_DELTA(force[0], expected_force_x, 1e-6);
            TS_ASSERT_DELTA(force[1], 0.0, 1e-10);
        }
    }

    /**
     * Test that the Morse and linear forces agree for small deviations
     * from equilibrium, and diverge for large deviations.
     */
    void TestSemMorseAndLinearForceAgreement()
    {
        SemForce<2> morse_force;
        SemLinearForce<2> linear_force;

        double u0 = 1.0;
        double rho = 5.0;
        double r_eq = 0.2;

        // For very small deviations about equilibrium, the Morse force should
        // agree with the tangent harmonic spring.
        double kappa = 8.0 * rho * rho * u0 / (r_eq * r_eq);
        for (double delta = -1e-4; delta <= 1e-4; delta += 5e-5)
        {
            double r = r_eq + delta;
            if (r <= 0.0) continue;

            c_vector<double, 2> vec_a_to_b;
            vec_a_to_b[0] = r;
            vec_a_to_b[1] = 0.0;
            double dist_sq = r * r;

            c_vector<double, 2> morse_f = morse_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            c_vector<double, 2> linear_f = linear_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            double expected_tangent_force = kappa * delta;

            TS_ASSERT_DELTA(linear_f[0], expected_tangent_force, 1e-10);
            TS_ASSERT_DELTA(morse_f[0], expected_tangent_force, 1e-2);
            TS_ASSERT_DELTA(linear_f[1], 0.0, 1e-10);
            TS_ASSERT_DELTA(morse_f[1], 0.0, 1e-10);
        }

        // At large compression (r = 0.1, 50% of r_eq), forces should diverge.
        // The Morse exponential should give much stronger repulsion than the linear spring.
        {
            double r = 0.1;
            c_vector<double, 2> vec_a_to_b;
            vec_a_to_b[0] = r;
            vec_a_to_b[1] = 0.0;
            double dist_sq = r * r;

            c_vector<double, 2> morse_f = morse_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);
            c_vector<double, 2> linear_f = linear_force.CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, r_eq);

            // Both should be repulsive (negative x)
            TS_ASSERT_LESS_THAN(morse_f[0], 0.0);
            TS_ASSERT_LESS_THAN(linear_f[0], 0.0);

            // Morse should be more strongly repulsive than linear at large compression
            TS_ASSERT_LESS_THAN(morse_f[0], linear_f[0]);
        }
    }

    c_vector<double, 3> CalculateExpectedSemMorseForce3d(const c_vector<double, 3>& rVectorAtoB,
                                                         double u0,
                                                         double rho,
                                                         double rEq)
    {
        const double distance_sq = inner_prod(rVectorAtoB, rVectorAtoB);
        const double r_eq_sq = rEq * rEq;
        const double s = distance_sq / r_eq_sq;
        const double exp_rho = std::exp(rho * (1.0 - s));
        const double coefficient = (4.0 * rho * u0 / r_eq_sq) * (exp_rho - exp_rho * exp_rho);

        return coefficient * rVectorAtoB;
    }

    /**
     * Test that CalculateForceBetweenNodes uses the intra-cellular parameter set
     * for two 3D SEM nodes in the same element.
     */
    void TestSemForceCalculateForceBetweenNodesUsesIntraParametersIn3d()
    {
        EXIT_IF_PARALLEL; // SEM is not parallel-ready
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 0.18, 0.24, 0.0));
        nodes[0]->AddElement(0u);
        nodes[1]->AddElement(0u);

        std::vector<Node<3>*> element_nodes;
        element_nodes.push_back(nodes[0]);
        element_nodes.push_back(nodes[1]);

        std::vector<SemElement<3>*> elements;
        elements.push_back(new SemElement<3>(0u, element_nodes));

        SemMesh<3> mesh(nodes, elements);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements());
        SemBasedCellPopulation<3> cell_population(mesh, cells);

        SemForce<3> force;
        force.SetIntraWellDepth(2.0);
        force.SetIntraScalingFactor(3.0);
        force.SetIntraEquilibriumDistance(0.2);
        force.SetIntraCutOffDistance(1.0);
        force.SetInterWellDepth(7.0);
        force.SetInterScalingFactor(11.0);
        force.SetInterEquilibriumDistance(0.8);
        force.SetInterCutOffDistance(1.0);

        const c_vector<double, 3> vec_a_to_b = nodes[1]->rGetLocation() - nodes[0]->rGetLocation();
        const c_vector<double, 3> expected_force = CalculateExpectedSemMorseForce3d(vec_a_to_b, 2.0, 3.0, 0.2);
        const c_vector<double, 3> force_between_nodes = force.CalculateForceBetweenNodes(0u, 1u, cell_population);

        TS_ASSERT_DELTA(force_between_nodes[0], expected_force[0], 1e-10);
        TS_ASSERT_DELTA(force_between_nodes[1], expected_force[1], 1e-10);
        TS_ASSERT_DELTA(force_between_nodes[2], expected_force[2], 1e-10);
    }

    /**
     * Test that CalculateForceBetweenNodes uses the inter-cellular parameter set
     * for two 3D SEM nodes in different elements.
     */
    void TestSemForceCalculateForceBetweenNodesUsesInterParametersIn3d()
    {
        EXIT_IF_PARALLEL; // SEM is not parallel-ready
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 0.18, 0.24, 0.0));
        nodes[0]->AddElement(0u);
        nodes[1]->AddElement(1u);

        std::vector<Node<3>*> element_0_nodes;
        element_0_nodes.push_back(nodes[0]);

        std::vector<Node<3>*> element_1_nodes;
        element_1_nodes.push_back(nodes[1]);

        std::vector<SemElement<3>*> elements;
        elements.push_back(new SemElement<3>(0u, element_0_nodes));
        elements.push_back(new SemElement<3>(1u, element_1_nodes));

        SemMesh<3> mesh(nodes, elements);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements());
        SemBasedCellPopulation<3> cell_population(mesh, cells);

        SemForce<3> force;
        force.SetIntraWellDepth(2.0);
        force.SetIntraScalingFactor(3.0);
        force.SetIntraEquilibriumDistance(0.2);
        force.SetIntraCutOffDistance(1.0);
        force.SetInterWellDepth(1.5);
        force.SetInterScalingFactor(4.0);
        force.SetInterEquilibriumDistance(0.25);
        force.SetInterCutOffDistance(1.0);

        const c_vector<double, 3> vec_a_to_b = nodes[1]->rGetLocation() - nodes[0]->rGetLocation();
        const c_vector<double, 3> expected_force = CalculateExpectedSemMorseForce3d(vec_a_to_b, 1.5, 4.0, 0.25);
        const c_vector<double, 3> force_between_nodes = force.CalculateForceBetweenNodes(0u, 1u, cell_population);

        TS_ASSERT_DELTA(force_between_nodes[0], expected_force[0], 1e-10);
        TS_ASSERT_DELTA(force_between_nodes[1], expected_force[1], 1e-10);
        TS_ASSERT_DELTA(force_between_nodes[2], expected_force[2], 1e-10);
    }

    /**
     * Test cutoff behaviour in 3D, including below/at/above cutoff and separate
     * intra- and inter-cellular cutoffs.
     */
    void TestSemForceCalculateForceBetweenNodesCutoffsIn3d()
    {
        EXIT_IF_PARALLEL; // SEM is not parallel-ready
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 0.18, 0.24, 0.0));
        nodes[0]->AddElement(0u);
        nodes[1]->AddElement(1u);

        std::vector<Node<3>*> element_0_nodes;
        element_0_nodes.push_back(nodes[0]);

        std::vector<Node<3>*> element_1_nodes;
        element_1_nodes.push_back(nodes[1]);

        std::vector<SemElement<3>*> elements;
        elements.push_back(new SemElement<3>(0u, element_0_nodes));
        elements.push_back(new SemElement<3>(1u, element_1_nodes));

        SemMesh<3> mesh(nodes, elements);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements());
        SemBasedCellPopulation<3> cell_population(mesh, cells);

        SemForce<3> force;
        force.SetInterWellDepth(1.0);
        force.SetInterScalingFactor(5.0);
        force.SetInterEquilibriumDistance(0.2);
        force.SetInterCutOffDistance(0.31);

        c_vector<double, 3> force_below_cutoff = force.CalculateForceBetweenNodes(0u, 1u, cell_population);
        TS_ASSERT_LESS_THAN(0.0, norm_2(force_below_cutoff));

        nodes[1]->rGetModifiableLocation()[0] = 0.3;
        nodes[1]->rGetModifiableLocation()[1] = 0.0;
        c_vector<double, 3> force_at_cutoff;
        force.SetInterCutOffDistance(0.3);
        force_at_cutoff = force.CalculateForceBetweenNodes(0u, 1u, cell_population);
        TS_ASSERT_DELTA(norm_2(force_at_cutoff), 0.0, 1e-12);

        nodes[1]->rGetModifiableLocation()[0] = 0.31;
        c_vector<double, 3> force_above_cutoff = force.CalculateForceBetweenNodes(0u, 1u, cell_population);
        TS_ASSERT_DELTA(norm_2(force_above_cutoff), 0.0, 1e-12);

        nodes[1]->rGetModifiableLocation()[0] = 0.18;
        nodes[1]->rGetModifiableLocation()[1] = 0.24;
        force.SetInterCutOffDistance(0.31);
        c_vector<double, 3> inter_force = force.CalculateForceBetweenNodes(0u, 1u, cell_population);
        TS_ASSERT_LESS_THAN(0.0, norm_2(inter_force));

        nodes[1]->rGetContainingElementIndices().clear();
        nodes[1]->AddElement(0u);
        force.SetIntraCutOffDistance(0.29);
        c_vector<double, 3> intra_force = force.CalculateForceBetweenNodes(0u, 1u, cell_population);
        TS_ASSERT_DELTA(norm_2(intra_force), 0.0, 1e-12);
    }

    /**
     * Test SemForce and SemLinearForce with a real SemBasedCellPopulation.
     * Uses a single-element mesh and verifies forces are finite.
     */
    void TestSemForceWithPopulation()
    {
        EXIT_IF_PARALLEL; // SEM is not parallel-ready
        // Create a single-element SEM mesh with nodes on a 5x5 grid, diameter 0.5
        SemSingleElementMeshGenerator<2> generator({5, 5}, 0.5);
        auto p_mesh = generator.GetMesh();

        c_vector<double, 4> domain;
        domain[0] = -1.0;
        domain[1] =  2.0;
        domain[2] = -1.0;
        domain[3] =  2.0;
        p_mesh->SetUpBoxCollection(0.5, domain);

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 25);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.0);

        // Populate node pairs via box collection
        cell_population.Update(false);

        // Clear forces
        for (unsigned i = 0; i < cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Apply the Morse force
        SemForce<2> morse_force;
        morse_force.AddForceContribution(cell_population);

        // Check that Morse forces are finite (no NaN or Inf)
        for (unsigned i = 0; i < cell_population.GetNumNodes(); i++)
        {
            for (unsigned d = 0; d < 2; d++)
            {
                TS_ASSERT(!std::isnan(cell_population.GetNode(i)->rGetAppliedForce()[d]));
                TS_ASSERT(!std::isinf(cell_population.GetNode(i)->rGetAppliedForce()[d]));
            }
        }

        // Clear Morse forces before applying the linear force
        for (unsigned i = 0; i < cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Apply the linear force with the same (default) parameters
        SemLinearForce<2> linear_force;
        linear_force.AddForceContribution(cell_population);

        // Check that linear forces are also finite
        for (unsigned i = 0; i < cell_population.GetNumNodes(); i++)
        {
            for (unsigned d = 0; d < 2; d++)
            {
                TS_ASSERT(!std::isnan(cell_population.GetNode(i)->rGetAppliedForce()[d]));
                TS_ASSERT(!std::isinf(cell_population.GetNode(i)->rGetAppliedForce()[d]));
            }
        }
    }

    /**
     * Test SemForce with a multi-element population to verify inter-cellular
     * forces are computed (using different parameters from intra-cellular).
     */
    void TestSemForceWithMultiElementPopulation()
    {
        EXIT_IF_PARALLEL; // SEM is not parallel-ready
        // Create a 2-element SEM mesh: two cells side by side
        SemMultiElementMeshGenerator<2> generator({5, 5}, {2, 1}, 0.5);
        auto p_mesh = generator.GetMesh();

        c_vector<double, 4> domain;
        domain[0] = -1.0;
        domain[1] =  3.0;
        domain[2] = -1.0;
        domain[3] =  2.0;
        p_mesh->SetUpBoxCollection(0.5, domain);

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 50);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.0);

        // Populate node pairs
        cell_population.Update(false);

        // Clear forces
        for (unsigned i = 0; i < cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Apply force with distinct inter/intra parameters to ensure both code paths are exercised
        SemForce<2> force;
        force.SetIntraWellDepth(1.0);
        force.SetIntraEquilibriumDistance(0.15);
        force.SetInterWellDepth(0.5);
        force.SetInterEquilibriumDistance(0.25);
        force.AddForceContribution(cell_population);

        // Check forces are finite
        bool any_nonzero = false;
        for (unsigned i = 0; i < cell_population.GetNumNodes(); i++)
        {
            for (unsigned d = 0; d < 2; d++)
            {
                double f = cell_population.GetNode(i)->rGetAppliedForce()[d];
                TS_ASSERT(!std::isnan(f));
                TS_ASSERT(!std::isinf(f));
                if (std::abs(f) > 1e-10)
                {
                    any_nonzero = true;
                }
            }
        }
        // At least some forces should be non-zero
        TS_ASSERT(any_nonzero);
    }

    /**
     * Test that the SEM Gaussian random force uses the SEM Langevin scaling.
     */
    void TestSemGaussianRandomForceWithPopulation()
    {
        EXIT_IF_PARALLEL; // SEM is not parallel-ready
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        SemSingleElementMeshGenerator<2> generator({3, 3}, 0.5);
        auto p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);

        SemGaussianRandomForce<2> force;
        force.SetDiffusionConstant(0.25);
        TS_ASSERT_DELTA(force.GetDiffusionConstant(), 0.25, 1e-12);

        RandomNumberGenerator::Instance()->Reseed(123u);
        force.AddForceContribution(cell_population);

        std::vector<c_vector<double, 2> > forces_with_lower_diffusion(cell_population.GetNumNodes());
        for (unsigned node_index = 0; node_index < cell_population.GetNumNodes(); ++node_index)
        {
            forces_with_lower_diffusion[node_index] = cell_population.GetNode(node_index)->rGetAppliedForce();
            cell_population.GetNode(node_index)->ClearAppliedForce();
        }

        force.SetDiffusionConstant(1.0);
        RandomNumberGenerator::Instance()->Reseed(123u);
        force.AddForceContribution(cell_population);

        bool any_nonzero = false;
        for (unsigned node_index = 0; node_index < cell_population.GetNumNodes(); ++node_index)
        {
            const c_vector<double, 2>& r_force_with_higher_diffusion = cell_population.GetNode(node_index)->rGetAppliedForce();
            for (unsigned dim = 0; dim < 2; ++dim)
            {
                TS_ASSERT_DELTA(r_force_with_higher_diffusion[dim], 2.0 * forces_with_lower_diffusion[node_index][dim], 1e-12);
                if (std::abs(r_force_with_higher_diffusion[dim]) > 1e-12)
                {
                    any_nonzero = true;
                }
            }
        }
        TS_ASSERT(any_nonzero);
    }

    /**
     * The noise can be cooled towards zero over a window of simulation time, so that an anneal, a
     * ramped quench and a subsequent noise-free relaxation all happen within one call to Solve().
     */
    void TestSemRandomForceCoolingWindow()
    {
        EXIT_IF_PARALLEL; // SEM is not parallel-ready
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10.0, 10);

        SemGaussianRandomForce<2> force;
        force.SetDiffusionConstant(0.8);

        // With no cooling window set, the diffusion constant never changes
        TS_ASSERT_DELTA(force.GetCurrentDiffusionConstant(), 0.8, 1e-12);
        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_DELTA(force.GetCurrentDiffusionConstant(), 0.8, 1e-12);

        // Cool from t = 2 to t = 6: full strength before, linearly down across, nothing after
        force.SetCoolingWindow(2.0, 6.0);
        TS_ASSERT_DELTA(force.GetCurrentDiffusionConstant(), 0.8, 1e-12);

        SimulationTime::Instance()->IncrementTimeOneStep(); // t = 2
        TS_ASSERT_DELTA(force.GetCurrentDiffusionConstant(), 0.8, 1e-12);

        SimulationTime::Instance()->IncrementTimeOneStep(); // t = 3, a quarter of the way through
        TS_ASSERT_DELTA(force.GetCurrentDiffusionConstant(), 0.6, 1e-12);

        SimulationTime::Instance()->IncrementTimeOneStep(); // t = 4, halfway
        TS_ASSERT_DELTA(force.GetCurrentDiffusionConstant(), 0.4, 1e-12);

        while (SimulationTime::Instance()->GetTime() < 6.0)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
        }
        TS_ASSERT_DELTA(force.GetCurrentDiffusionConstant(), 0.0, 1e-12);

        SimulationTime::Instance()->IncrementTimeOneStep(); // t = 7, past the end of the window
        TS_ASSERT_DELTA(force.GetCurrentDiffusionConstant(), 0.0, 1e-12);

        // The value set by SetDiffusionConstant() is left alone throughout
        TS_ASSERT_DELTA(force.GetDiffusionConstant(), 0.8, 1e-12);

        TS_ASSERT_THROWS_THIS(force.SetCoolingWindow(6.0, 2.0),
            "AbstractSemRandomForce: the cooling window must not end before it starts");
    }

    /**
     * A cooled force must apply correspondingly smaller forces, and none at all once the window
     * has passed.
     */
    void TestSemRandomForceCoolingWindowScalesTheAppliedForce()
    {
        EXIT_IF_PARALLEL; // SEM is not parallel-ready
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(8.0, 8);

        SemSingleElementMeshGenerator<2> generator({3, 3}, 0.5);
        auto p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);

        SemGaussianRandomForce<2> force;
        force.SetDiffusionConstant(1.0);
        force.SetCoolingWindow(0.0, 4.0);

        // At t = 1 the window is a quarter spent, so D is 0.75 of its full value and the force,
        // which goes as the square root of D, is scaled by sqrt(0.75)
        SimulationTime::Instance()->IncrementTimeOneStep();
        RandomNumberGenerator::Instance()->Reseed(123u);
        force.AddForceContribution(cell_population);

        std::vector<c_vector<double, 2> > cooled_forces(cell_population.GetNumNodes());
        for (unsigned i = 0; i < cell_population.GetNumNodes(); ++i)
        {
            cooled_forces[i] = cell_population.GetNode(i)->rGetAppliedForce();
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        SemGaussianRandomForce<2> uncooled_force;
        uncooled_force.SetDiffusionConstant(1.0);
        RandomNumberGenerator::Instance()->Reseed(123u);
        uncooled_force.AddForceContribution(cell_population);

        for (unsigned i = 0; i < cell_population.GetNumNodes(); ++i)
        {
            const c_vector<double, 2>& r_full = cell_population.GetNode(i)->rGetAppliedForce();
            for (unsigned dim = 0; dim < 2; ++dim)
            {
                TS_ASSERT_DELTA(cooled_forces[i][dim], sqrt(0.75) * r_full[dim], 1e-12);
            }
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Once the window has passed there is no noise at all
        while (SimulationTime::Instance()->GetTime() < 5.0)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
        }
        force.AddForceContribution(cell_population);
        for (unsigned i = 0; i < cell_population.GetNumNodes(); ++i)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), 0.0, 1e-12);
        }
    }

    /**
     * Test the SEM spatially correlated random force with a real SEM population.
     */
    void TestSemSpatiallyCorrelatedRandomForceWithPopulation()
    {
        EXIT_IF_PARALLEL; // SEM is not parallel-ready
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        SemSingleElementMeshGenerator<2> generator({3, 3}, 0.5);
        auto p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.0);

        SemSpatiallyCorrelatedRandomForce<2> force;
        force.SetDiffusionConstant(0.25);
        force.SetCorrelationLength(0.5);
        force.SetLowerCorner({{-1.0, -1.0}});
        force.SetUpperCorner({{1.0, 1.0}});
        force.SetPeriodicity({{false, false}});
        force.SetRandomSeed(123u);

        TS_ASSERT_DELTA(force.GetDiffusionConstant(), 0.25, 1e-12);
        TS_ASSERT_DELTA(force.GetCorrelationLength(), 0.5, 1e-12);
        TS_ASSERT_EQUALS(force.GetLowerCorner()[0], -1.0);
        TS_ASSERT_EQUALS(force.GetUpperCorner()[1], 1.0);
        TS_ASSERT_EQUALS(force.GetPeriodicity()[0], false);

        RandomNumberGenerator::Instance()->Reseed(123u);
        force.AddForceContribution(cell_population);

        bool any_nonzero = false;
        for (unsigned node_index = 0; node_index < cell_population.GetNumNodes(); ++node_index)
        {
            for (unsigned dim = 0; dim < 2; ++dim)
            {
                const double force_component = cell_population.GetNode(node_index)->rGetAppliedForce()[dim];
                TS_ASSERT(!std::isnan(force_component));
                TS_ASSERT(!std::isinf(force_component));
                if (std::abs(force_component) > 1e-12)
                {
                    any_nonzero = true;
                }
            }
        }
        TS_ASSERT(any_nonzero);
    }

    /**
     * Test that SemForce throws an exception when used with a non-SEM population.
     */
    void TestSemForceWithWrongPopulationType()
    {
        // Create a simple NodeBasedCellPopulation
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.1, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.0);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        SemForce<2> sem_force;
        TS_ASSERT_THROWS_THIS(sem_force.AddForceContribution(cell_population),
                "SemForce is to be used with a SemBasedCellPopulation only");

        SemLinearForce<2> sem_linear_force;
        TS_ASSERT_THROWS_THIS(sem_linear_force.AddForceContribution(cell_population),
                "SemLinearForce is to be used with a SemBasedCellPopulation only");

        SemGaussianRandomForce<2> sem_gaussian_random_force;
        TS_ASSERT_THROWS_THIS(sem_gaussian_random_force.AddForceContribution(cell_population),
                "AbstractSemRandomForce is to be used with a SemBasedCellPopulation only");

        SemSpatiallyCorrelatedRandomForce<2> sem_spatially_correlated_random_force;
        TS_ASSERT_THROWS_THIS(sem_spatially_correlated_random_force.AddForceContribution(cell_population),
                "AbstractSemRandomForce is to be used with a SemBasedCellPopulation only");

        for (unsigned i = 0; i < nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    /**
     * Test archiving of SemForce.
     */
    void TestSemForceArchiving()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SemForce.arch";

        {
            SemForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetIntraWellDepth(2.5);
            force.SetIntraScalingFactor(3.0);
            force.SetIntraEquilibriumDistance(0.15);
            force.SetIntraCutOffDistance(0.6);
            force.SetInterWellDepth(1.5);
            force.SetInterScalingFactor(4.0);
            force.SetInterEquilibriumDistance(0.25);
            force.SetInterCutOffDistance(0.7);

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
            TS_ASSERT_DELTA((static_cast<SemForce<2>*>(p_force))->GetIntraWellDepth(), 2.5, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemForce<2>*>(p_force))->GetIntraScalingFactor(), 3.0, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemForce<2>*>(p_force))->GetIntraEquilibriumDistance(), 0.15, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemForce<2>*>(p_force))->GetIntraCutOffDistance(), 0.6, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemForce<2>*>(p_force))->GetInterWellDepth(), 1.5, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemForce<2>*>(p_force))->GetInterScalingFactor(), 4.0, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemForce<2>*>(p_force))->GetInterEquilibriumDistance(), 0.25, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemForce<2>*>(p_force))->GetInterCutOffDistance(), 0.7, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    /**
     * Test archiving of SemLinearForce.
     */
    void TestSemLinearForceArchiving()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SemLinearForce.arch";

        {
            SemLinearForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetIntraWellDepth(2.5);
            force.SetIntraScalingFactor(3.0);
            force.SetIntraEquilibriumDistance(0.15);
            force.SetIntraCutOffDistance(0.6);

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
            TS_ASSERT_DELTA((static_cast<SemLinearForce<2>*>(p_force))->GetIntraWellDepth(), 2.5, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemLinearForce<2>*>(p_force))->GetIntraScalingFactor(), 3.0, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemLinearForce<2>*>(p_force))->GetIntraEquilibriumDistance(), 0.15, 1e-6);
            TS_ASSERT_DELTA((static_cast<SemLinearForce<2>*>(p_force))->GetIntraCutOffDistance(), 0.6, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    /**
     * Test archiving of SEM random forces.
     */
    void TestSemRandomForceArchiving()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler handler("archive", false);

        {
            const std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SemGaussianRandomForce.arch";

            {
                SemGaussianRandomForce<2> force;
                force.SetDiffusionConstant(0.25);

                std::ofstream ofs(archive_filename.c_str());
                boost::archive::text_oarchive output_arch(ofs);

                AbstractForce<2>* const p_force = &force;
                output_arch << p_force;
            }

            {
                AbstractForce<2>* p_force;

                std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
                boost::archive::text_iarchive input_arch(ifs);

                input_arch >> p_force;

                TS_ASSERT_DELTA((static_cast<SemGaussianRandomForce<2>*>(p_force))->GetDiffusionConstant(), 0.25, 1e-6);

                delete p_force;
            }
        }

        {
            const std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SemSpatiallyCorrelatedRandomForce.arch";

            {
                SemSpatiallyCorrelatedRandomForce<2> force;
                force.SetDiffusionConstant(0.5);
                force.SetCorrelationLength(0.25);
                force.SetLowerCorner({{-1.0, -2.0}});
                force.SetUpperCorner({{3.0, 4.0}});
                force.SetPeriodicity({{true, false}});
                force.SetSmallCorrelationLengthWarningThreshold(1e-9);

                std::ofstream ofs(archive_filename.c_str());
                boost::archive::text_oarchive output_arch(ofs);

                AbstractForce<2>* const p_force = &force;
                output_arch << p_force;
            }

            {
                AbstractForce<2>* p_force;

                std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
                boost::archive::text_iarchive input_arch(ifs);

                input_arch >> p_force;

                auto p_random_force = static_cast<SemSpatiallyCorrelatedRandomForce<2>*>(p_force);
                TS_ASSERT_DELTA(p_random_force->GetDiffusionConstant(), 0.5, 1e-6);
                TS_ASSERT_DELTA(p_random_force->GetCorrelationLength(), 0.25, 1e-6);
                TS_ASSERT_DELTA(p_random_force->GetLowerCorner()[0], -1.0, 1e-6);
                TS_ASSERT_DELTA(p_random_force->GetLowerCorner()[1], -2.0, 1e-6);
                TS_ASSERT_DELTA(p_random_force->GetUpperCorner()[0], 3.0, 1e-6);
                TS_ASSERT_DELTA(p_random_force->GetUpperCorner()[1], 4.0, 1e-6);
                TS_ASSERT_EQUALS(p_random_force->GetPeriodicity()[0], true);
                TS_ASSERT_EQUALS(p_random_force->GetPeriodicity()[1], false);
                TS_ASSERT_DELTA(p_random_force->GetSmallCorrelationLengthWarningThreshold(), 1e-9, 1e-12);

                delete p_force;
            }
        }
    }
};

#endif /*TESTFORCES_HPP_*/
