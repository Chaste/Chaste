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

#ifndef TESTDISCRETESYSTEMFORCECALCULATOR_HPP_
#define TESTDISCRETESYSTEMFORCECALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "OffLatticeSimulation.hpp"
#include "DiscreteSystemForceCalculator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "NumericFileComparison.hpp"
#include "FakePetscSetup.hpp"

class TestDiscreteSystemForceCalculator : public AbstractCellBasedTestSuite
{
public:

    void TestPrivateMethods()
    {
        // Set up a cell population
        HoneycombMeshGenerator mesh_generator(7, 5, 0, 2.0);
        MutableMesh<2,2>* p_mesh = mesh_generator.GetMesh();
        std::vector<unsigned> location_indices = mesh_generator.GetCellLocationIndices();

        CellsGenerator<FixedG1GenerationalCellCycleModel,2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create the force law and pass in to a std::list
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        std::vector<boost::shared_ptr<AbstractTwoBodyInteractionForce<2> > > force_collection;
        force_collection.push_back(p_force);

        // Create a force calculator
        DiscreteSystemForceCalculator calculator(cell_population, force_collection);

        unsigned node_index = 8;

        // Test GetNeighbouringNodeIndices
        std::set<unsigned> expected_node_indices;
        expected_node_indices.insert(1u);
        expected_node_indices.insert(2u);
        expected_node_indices.insert(7u);
        expected_node_indices.insert(9u);
        expected_node_indices.insert(15u);
        expected_node_indices.insert(16u);

        std::set<unsigned> neighbouring_node_indices = cell_population.GetNeighbouringNodeIndices(node_index);

        TS_ASSERT(neighbouring_node_indices == expected_node_indices);

        // Test CalculateFtAndFn
        double spring_stiffness = p_force->GetMeinekeSpringStiffness();
        double expected_ft = spring_stiffness*(cos(M_PI/12.0) + cos(5.0*M_PI/12.0) - cos(3.0*M_PI/12.0));
        double expected_fn = spring_stiffness*(sin(M_PI/12.0) + sin(5.0*M_PI/12.0) + sin(3.0*M_PI/12.0));
        std::vector<double> Ft_and_Fn = calculator.CalculateFtAndFn(node_index,M_PI/4.0);
        TS_ASSERT_DELTA(Ft_and_Fn[0], expected_ft, 1e-4);
        TS_ASSERT_DELTA(Ft_and_Fn[1], expected_fn, 1e-4);

        // Test GetSamplingAngles
        std::vector<double> expected_sampling_angles;
        double epsilon = calculator.mEpsilon;

        expected_sampling_angles.push_back(-M_PI +epsilon);
        expected_sampling_angles.push_back(-M_PI +epsilon);

        for (unsigned i=1; i<6; i++)
        {
            expected_sampling_angles.push_back(-M_PI +((double) i)*M_PI/3.0 - epsilon);
            expected_sampling_angles.push_back(-M_PI +((double) i)*M_PI/3.0 - epsilon);
            expected_sampling_angles.push_back(-M_PI +((double) i)*M_PI/3.0 + epsilon);
            expected_sampling_angles.push_back(-M_PI +((double) i)*M_PI/3.0 + epsilon);
        }
        expected_sampling_angles.push_back(M_PI - epsilon);
        expected_sampling_angles.push_back(M_PI - epsilon);

        std::vector<double> sampling_angles = calculator.GetSamplingAngles(node_index);
        for (unsigned i=0; i<sampling_angles.size(); i++)
        {
            // the sampling angles lie in the range (pi,pi]
            if (expected_sampling_angles[i] > M_PI)
            {
                expected_sampling_angles[i] -= 2*M_PI;
            }

            TS_ASSERT_DELTA(sampling_angles[i], expected_sampling_angles[i], 1e-6);
        }

        // Test GetLocalExtremum
        double expected_extremal_angle = -M_PI +M_PI/6.0;
        double calculated_extremal_angle = calculator.GetLocalExtremum(node_index, sampling_angles[1], sampling_angles[2]);

        TS_ASSERT_DELTA(calculated_extremal_angle, expected_extremal_angle, 1e-4);

        // Test GetExtremalAngles
        std::vector<double> calculated_extremal_angles = calculator.GetExtremalAngles(node_index, sampling_angles);

        // the extremal angles lie in the range (pi,pi]
        TS_ASSERT_DELTA(-M_PI + M_PI/6.0, calculated_extremal_angles[0], 1e-4);
        TS_ASSERT_DELTA(-M_PI + M_PI/3.0, calculated_extremal_angles[1], 1e-4);
        TS_ASSERT_DELTA(-M_PI + M_PI/2.0, calculated_extremal_angles[2], 1e-4);
        TS_ASSERT_DELTA(-M_PI + 2.0*M_PI/3.0, calculated_extremal_angles[3], 1e-4);
        TS_ASSERT_DELTA(-M_PI + 5.0*M_PI/6.0, calculated_extremal_angles[4], 1e-4);
    }

    void TestCalculateExtremalNormalForces()
    {
        // Set up a cell population
        HoneycombMeshGenerator mesh_generator(7, 5, 0, 2.0);
        MutableMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create the force law and pass in to a std::list
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        std::vector<boost::shared_ptr<AbstractTwoBodyInteractionForce<2> > > force_collection;
        force_collection.push_back(p_force);

        // Create a force calculator
        DiscreteSystemForceCalculator calculator(cell_population, force_collection);

        // Test CalculateExtremalNormalForces
        std::vector< std::vector<double> > calculated_results = calculator.CalculateExtremalNormalForces();
        TS_ASSERT_EQUALS(calculated_results.size(), 2u);
        TS_ASSERT_EQUALS(calculated_results[0].size(), p_mesh->GetNumNodes());
        TS_ASSERT_EQUALS(calculated_results[1].size(), p_mesh->GetNumNodes());

        double spring_stiffness = p_force->GetMeinekeSpringStiffness();
        double expected_minimum_interior = spring_stiffness*( 2.0*sin(M_PI/3.0) );
        double expected_maximum_interior = spring_stiffness*( sin(M_PI/6.0) + sin(M_PI/2.0) + sin(5.0*M_PI/6.0) );

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            if (!(p_mesh->GetNode(i)->IsBoundaryNode()))
            {
                TS_ASSERT_DELTA(calculated_results[0][i], expected_minimum_interior, 1e-4);
                TS_ASSERT_DELTA(calculated_results[1][i], expected_maximum_interior, 1e-4);
            }
            else
            {
                // Separate cases for the boundary nodes...
                if (i==29 || i==30 || i==31 || i==32 || i==33)
                {
                    TS_ASSERT_DELTA(calculated_results[0][i], 0.0, 1e-4);
                    TS_ASSERT_DELTA(calculated_results[1][i], expected_maximum_interior, 1e-4);
                }
                if (i==7 || i==20 || i==21)
                {
                    TS_ASSERT_DELTA(calculated_results[0][i], spring_stiffness*sin(M_PI/3.0), 1e-4);
                    TS_ASSERT_DELTA(calculated_results[1][i], expected_maximum_interior, 1e-4);

                }
                if (i==1 || i==2 || i==3 || i==4 || i==5)
                {
                    TS_ASSERT_DELTA(calculated_results[0][i], 0.0, 1e-4);
                    TS_ASSERT_DELTA(calculated_results[1][i], expected_maximum_interior, 1e-4);

                }
                if (i==13 || i==14 || i==27)
                {
                    TS_ASSERT_DELTA(calculated_results[0][i], expected_maximum_interior, 1e-4);
                    TS_ASSERT_DELTA(calculated_results[1][i], expected_maximum_interior, 1e-4);
                }
            }
        }
    }

    void TestCalculateWriteResultsToFile()
    {
        std::string output_directory = "TestDiscreteSystemForceCalculator";

        // Set up a cell population

        HoneycombMeshGenerator mesh_generator(7, 5, 0, 2.0);
        MutableMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create the force law and pass in to a std::list
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        std::vector<boost::shared_ptr<AbstractTwoBodyInteractionForce<2> > > force_collection;
        force_collection.push_back(p_force);

        // Create a force calculator
        DiscreteSystemForceCalculator calculator(cell_population, force_collection);

        // Test WriteResultsToFile
        calculator.WriteResultsToFile(output_directory);

        // Compare output with saved files of what they should look like
        OutputFileHandler handler(output_directory, false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizstress";

        NumericFileComparison node_velocities(results_file, "cell_based/test/data/TestDiscreteSystemForceCalculator/results.vizstress");
        TS_ASSERT(node_velocities.CompareFiles(1e-4));

        // Run a simulation to generate some results.viz<other things> files
        // so the visualizer can display the results.vizstress file.
        // (These lines are not actually necessary for generating results.vizstress)
        OffLatticeSimulation<2> simulator(cell_population);

        simulator.AddForce(p_force);

        simulator.SetEndTime(0.05);
        simulator.SetOutputDirectory(output_directory+"_rerun");
        simulator.Solve();
    }

    void TestAllCases()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA(mesh.GetAngleBetweenNodes(2,0), -0.75*M_PI, 1e-12);
        TS_ASSERT_DELTA(mesh.GetAngleBetweenNodes(2,1), -0.5*M_PI, 1e-12);

        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        std::vector<boost::shared_ptr<AbstractTwoBodyInteractionForce<2> > > force_collection;
        force_collection.push_back(p_force);

        DiscreteSystemForceCalculator calculator(cell_population, force_collection);

        double epsilon = 0.5*M_PI;
        calculator.SetEpsilon(epsilon);

        TS_ASSERT_THROWS_NOTHING(calculator.GetSamplingAngles(2));
    }
};

#endif /*TESTDISCRETESYSTEMFORCECALCULATOR_HPP_*/
