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

#ifndef TESTBUSKEFORCES_HPP_
#define TESTBUSKEFORCES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "BuskeAdhesiveForce.hpp"
#include "BuskeElasticForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FileComparison.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestBuskeForces : public AbstractCellBasedTestSuite
{
public:

    void TestBuskeAdhesiveForceMethods()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 10.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 20.0);
        mesh.GetNode(0)->SetRadius(1.0);
        mesh.GetNode(1)->SetRadius(1.0);
        mesh.GetNode(2)->SetRadius(1.0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Update();

        // Create force
        BuskeAdhesiveForce<2> buske_adhesive_force;

        // Test set/get methods
        TS_ASSERT_EQUALS(buske_adhesive_force.GetUseCutOffLength(), false);
        TS_ASSERT_DELTA(buske_adhesive_force.GetCutOffLength(), DBL_MAX, 1e-6);
        TS_ASSERT_DELTA(buske_adhesive_force.GetAdhesionEnergyParameter(), 0.2, 1e-6);

        buske_adhesive_force.SetCutOffLength(1.5);
        buske_adhesive_force.SetAdhesionEnergyParameter(1.0);

        TS_ASSERT_EQUALS(buske_adhesive_force.GetUseCutOffLength(), true);
        TS_ASSERT_DELTA(buske_adhesive_force.GetCutOffLength(), 1.5, 1e-6);
        TS_ASSERT_DELTA(buske_adhesive_force.GetAdhesionEnergyParameter(), 1.0, 1e-6);

        buske_adhesive_force.SetCutOffLength(5.0);
        buske_adhesive_force.SetAdhesionEnergyParameter(200);

        // Test node force calculation
        for (unsigned step=0; step<40; step++)
        {
            // Move nodes close together
            double separation = 4.0 - (double)step/10.0;
            cell_population.GetNode(1)->rGetModifiableLocation()[0] = separation;

            // Reset the vector of node forces
            for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
            {
                cell_population.GetNode(i)->ClearAppliedForce();
            }

            buske_adhesive_force.AddForceContribution(cell_population);

            // Test forces on nodes
            double analytical_force_magnitude = 0.0;
            if (separation < 2.0)
            {
                analytical_force_magnitude = 100.0*M_PI*separation;
            }

            TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], analytical_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);

            TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], -analytical_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);

            TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0, 1e-4);
        }

        // Clean up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestBuskeElasticForceMethods()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 10.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 20.0);
        mesh.GetNode(0)->SetRadius(1.0);
        mesh.GetNode(1)->SetRadius(1.0);
        mesh.GetNode(2)->SetRadius(1.0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Update();

        // Create force
        BuskeElasticForce<2> buske_elastic_force;

        // Test set/get methods
        TS_ASSERT_EQUALS(buske_elastic_force.GetUseCutOffLength(), false);
        TS_ASSERT_DELTA(buske_elastic_force.GetCutOffLength(), DBL_MAX, 1e-6);
        TS_ASSERT_DELTA(buske_elastic_force.GetDeformationEnergyParameter(), 4.0/(3.0*5.0), 1e-6);

        buske_elastic_force.SetCutOffLength(1.5);
        buske_elastic_force.SetDeformationEnergyParameter(1.0);

        TS_ASSERT_EQUALS(buske_elastic_force.GetUseCutOffLength(), true);
        TS_ASSERT_DELTA(buske_elastic_force.GetCutOffLength(), 1.5, 1e-6);
        TS_ASSERT_DELTA(buske_elastic_force.GetDeformationEnergyParameter(), 1.0, 1e-6);

        buske_elastic_force.SetCutOffLength(5.0);
        buske_elastic_force.SetDeformationEnergyParameter(4.0/3.0);

        // Test node force calculation
        for (unsigned step=0; step<40; step++)
        {
            // Move nodes close together
            double separation = 4.0 - (double)step/10.0;
            cell_population.GetNode(1)->rGetModifiableLocation()[0] = separation;

            // Reset the vector of node forces
            for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
            {
                cell_population.GetNode(i)->ClearAppliedForce();
            }

            buske_elastic_force.AddForceContribution(cell_population);

            // Test forces on nodes
            double analytical_force_magnitude = 0.0;//100.0*M_PI*separation;
            if (separation < 2.0)
            {
                analytical_force_magnitude -= 3.0/4.0*pow((2.0-separation),1.5)*sqrt(0.5);
            }

            ///\todo test force calculation (#1764)
            TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], analytical_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);

            TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], -analytical_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);

            TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0, 1e-4);
        }

        // Clean up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestBuskeMixedForceMethods()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        HoneycombMeshGenerator generator(2, 1, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 2.5);
        mesh.GetNode(0)->SetRadius(1.0);
        mesh.GetNode(1)->SetRadius(1.0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Update();

        // Create force
        BuskeAdhesiveForce<2> buske_adhesive_force;
        BuskeElasticForce<2> buske_elastic_force;
        BuskeCompressionForce<2> buske_compression_force;

        // Test node force calculation

        // Initialise a vector of node forces
        std::vector<c_vector<double, 2> > node_adhesive_force;
        std::vector<c_vector<double, 2> > node_elastic_force;
        std::vector<c_vector<double, 2> > node_compression_force;

        node_adhesive_force.reserve(cell_population.GetNumNodes());
        node_elastic_force.reserve(cell_population.GetNumNodes());
        node_compression_force.reserve(cell_population.GetNumNodes());

        for (unsigned step=0; step<40; step++)
        {
            // Move nodes close together
            double separation = 4.0 - (double)step/10.0;
            cell_population.GetNode(1)->rGetModifiableLocation()[0] = separation;

            // Reset the vector of node forces
            for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
            {
                cell_population.GetNode(i)->ClearAppliedForce();
            }

            buske_adhesive_force.AddForceContribution(cell_population);
            buske_elastic_force.AddForceContribution(cell_population);
            buske_compression_force.AddForceContribution(cell_population);

            // Test forces
            ///\todo test something!
        }
    }

    void TestBuskeCompressionForceMethods()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        HoneycombMeshGenerator generator(2, 1, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 2.5);
        mesh.GetNode(0)->SetRadius(1.0);
        mesh.GetNode(1)->SetRadius(1.0);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Update();

        // Create force
        BuskeCompressionForce<2> buske_compression_force;

        // Test set/get methods
        TS_ASSERT_DELTA(buske_compression_force.GetCompressionEnergyParameter(), 5.0, 1e-6);
        buske_compression_force.SetCompressionEnergyParameter(15.0);

        TS_ASSERT_DELTA(buske_compression_force.GetCompressionEnergyParameter(), 15.0, 1e-6);
        buske_compression_force.SetCompressionEnergyParameter(1.0);

        // Test node force calculation
        for (unsigned step=0; step<40; step++)
        {
            // Move nodes close together
            double separation = 4.0 - (double)step/10.0;
            cell_population.GetNode(1)->rGetModifiableLocation()[0] = separation;

            // Reset the vector of node forces
            for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
            {
                cell_population.GetNode(i)->ClearAppliedForce();
            }

            buske_compression_force.AddForceContribution(cell_population);

            // Test forces on nodes
            double analytical_force_magnitude = 0.0;
            if (separation < 2.0)
            {
                analytical_force_magnitude = M_PI/6.0/5.0*(5.0-M_PI/3.0*(4.0-pow((1.0-separation/2.0),2.0)*(2.0-separation/2.0)))
                                                          *(3.0/4.0*pow(separation,2.0)-4.0*separation+5.0);
            }

            TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], -analytical_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);

            TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], analytical_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);
        }
    }

    void TestBuskeCompressionForceWithMultipleCells()
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create triangle mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true,  0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, -1.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 2.5);
        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->SetRadius(1.0);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Update();

        // Create force
        BuskeCompressionForce<2> buske_compression_force;

        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        buske_compression_force.AddForceContribution(cell_population);

        // Test forces on nodes

        // This node should only move in the x-direction
        if (PetscTools::AmMaster())    // All cells should start on the master.
        {
            unsigned zero_node_index = 0;
            unsigned first_node_index = PetscTools::GetNumProcs();
            unsigned second_node_index = 2 * PetscTools::GetNumProcs();
            TS_ASSERT_DELTA(cell_population.GetNode(zero_node_index)->rGetAppliedForce()[0], -0.6514, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(zero_node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(first_node_index)->rGetAppliedForce()[0], 0.2894, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(first_node_index)->rGetAppliedForce()[1], -0.2894, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(second_node_index)->rGetAppliedForce()[0], 0.2894, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(second_node_index)->rGetAppliedForce()[1], 0.2894, 1e-4);
        }

        // When the node-only mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestForceOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestNotForReleaseForcesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with BuskeAdhesiveForce
        BuskeAdhesiveForce<2> buske_adhesive_force;
        buske_adhesive_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(buske_adhesive_force.GetIdentifier(), "BuskeAdhesiveForce-2");

        out_stream buske_force_parameter_file = output_file_handler.OpenOutputFile("buske_results.parameters");
        buske_adhesive_force.OutputForceParameters(buske_force_parameter_file);

        // Test with BuskeElasticForce
        BuskeElasticForce<2> buske_elastic_force;
        buske_elastic_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(buske_elastic_force.GetIdentifier(), "BuskeElasticForce-2");

        buske_elastic_force.OutputForceParameters(buske_force_parameter_file);
        buske_force_parameter_file->close();

        {
            FileFinder generated = output_file_handler.FindFile("buske_results.parameters");
            FileFinder reference("cell_based/test/data/TestForcesNotForRelease/buske_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with BuskeCompressionForce
        BuskeCompressionForce<2> buske_compression_force;
        TS_ASSERT_EQUALS(buske_compression_force.GetIdentifier(), "BuskeCompressionForce-2");

        out_stream buske_compression_force_parameter_file = output_file_handler.OpenOutputFile("buske_compression_results.parameters");
        buske_compression_force.OutputForceParameters(buske_compression_force_parameter_file);
        buske_compression_force_parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code.
            FileFinder generated = output_file_handler.FindFile("buske_compression_results.parameters");
            FileFinder reference("cell_based/test/data/TestForcesNotForRelease/buske_compression_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestBuskeAdhesiveForceArchiving()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "BuskeAdhesiveForce.arch";

        {
            BuskeAdhesiveForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetCutOffLength(1.7);
            force.SetAdhesionEnergyParameter(12.0);

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
            TS_ASSERT_EQUALS(dynamic_cast<BuskeAdhesiveForce<2>*>(p_force)->GetUseCutOffLength(), true);
            TS_ASSERT_DELTA(dynamic_cast<BuskeAdhesiveForce<2>*>(p_force)->GetCutOffLength(), 1.7, 1e-6);
            TS_ASSERT_DELTA(dynamic_cast<BuskeAdhesiveForce<2>*>(p_force)->GetAdhesionEnergyParameter(), 12.0, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestBuskeElasticForceArchiving()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "BuskeElasticForce.arch";

        {
            BuskeElasticForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetCutOffLength(1.7);
            force.SetDeformationEnergyParameter(13.0);

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
            TS_ASSERT_EQUALS(dynamic_cast<BuskeElasticForce<2>*>(p_force)->GetUseCutOffLength(), true);
            TS_ASSERT_DELTA(dynamic_cast<BuskeElasticForce<2>*>(p_force)->GetCutOffLength(), 1.7, 1e-6);
            TS_ASSERT_DELTA(dynamic_cast<BuskeElasticForce<2>*>(p_force)->GetDeformationEnergyParameter(), 13.0, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestBuskeCompressionForceArchiving()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "BuskeCompressionForce.arch";

        {
            BuskeCompressionForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetCompressionEnergyParameter(14.0);

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
            TS_ASSERT_DELTA(dynamic_cast<BuskeCompressionForce<2>*>(p_force)->GetCompressionEnergyParameter(), 14.0, 1e-6);

            // Tidy up
            delete p_force;
        }
    }
};

#endif /*TESTBUSKEFORCES_HPP_*/
