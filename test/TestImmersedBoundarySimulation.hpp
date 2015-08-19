/*

Copyright (c) 2005-2014, University of Oxford.
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

// Needed for test framework
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// Needed for Immersed Boundary simulations
#include <fftw3.h>

// Includes from trunk
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "StochasticDurationCellCycleModel.hpp"


// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryMeshWriter.hpp"
#include "ImmersedBoundaryMeshReader.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "SuperellipseGenerator.hpp"

// Immersed boundary forces
#include "ImmersedBoundaryElasticityForce.hpp"

#include "Debug.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundarySimulation : public AbstractCellBasedTestSuite
{
public:

    void xTestImmersedBoundaryMeshArchiving() throw(Exception)
    {
        // Create a vector of nodes forming a rectangle in (0,1)x(0,1)
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(1, true, 0.2, 0.1));
        nodes.push_back(new Node<2>(2, true, 0.3, 0.1));
        nodes.push_back(new Node<2>(3, true, 0.4, 0.1));
        nodes.push_back(new Node<2>(4, true, 0.4, 0.2));
        nodes.push_back(new Node<2>(5, true, 0.4, 0.3));
        nodes.push_back(new Node<2>(6, true, 0.4, 0.4));
        nodes.push_back(new Node<2>(7, true, 0.4, 0.5));
        nodes.push_back(new Node<2>(8, true, 0.4, 0.6));
        nodes.push_back(new Node<2>(9, true, 0.4, 0.7));
        nodes.push_back(new Node<2>(10, true, 0.4, 0.8));
        nodes.push_back(new Node<2>(11, true, 0.3, 0.8));
        nodes.push_back(new Node<2>(12, true, 0.2, 0.8));
        nodes.push_back(new Node<2>(13, true, 0.1, 0.8));
        nodes.push_back(new Node<2>(14, true, 0.1, 0.7));
        nodes.push_back(new Node<2>(15, true, 0.1, 0.6));
        nodes.push_back(new Node<2>(16, true, 0.1, 0.5));
        nodes.push_back(new Node<2>(17, true, 0.1, 0.4));
        nodes.push_back(new Node<2>(18, true, 0.1, 0.3));
        nodes.push_back(new Node<2>(19, true, 0.1, 0.2));

        // Create a vector of immersed boundary elements and create an element with the nodes above
        std::vector<ImmersedBoundaryElement<2,2>*> ib_element;
        ib_element.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));

        // Create a mesh with the nodes and elements vectors
        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_element, 32, 64);

        // Modify the fluid grid in one place to help test mesh writer and reader
        std::vector<std::vector<double> >& modifiable_fluid_grid = p_mesh->rGetModifiableFluidVelocityGridX();
        modifiable_fluid_grid[5][3] = 12.0;

        // Write the state of the immersed boundary mesh to file
        ImmersedBoundaryMeshWriter<2,2> mesh_writer("IBMeshOneSquareElement", "ib_mesh_one_square_element");
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }

    void xTestImmersedBoundaryMeshReading() throw(Exception)
    {
        // Load immersed boundary mesh
        ImmersedBoundaryMeshReader<2,2> mesh_reader("projects/ImmersedBoundary/test/mesh/ib_mesh_one_square_element");

        // Construct the immersed boundary mesh from the mesh reader
        ImmersedBoundaryMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Test the point that was set to 12 before mesh was written
        TS_ASSERT_DELTA(mesh.rGetFluidVelocityGridX()[5][3], 12.0, 1e-6);
    }

    void xTestRetrieveElementProperties() throw(Exception)
    {
        // Load immersed boundary mesh
        ImmersedBoundaryMeshReader<2,2> mesh_reader("projects/ImmersedBoundary/test/mesh/ib_mesh_one_square_element");

        // Construct the immersed boundary mesh from the mesh reader
        ImmersedBoundaryMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        ImmersedBoundaryMesh<2,2>* p_mesh = &mesh;

        /**
         * Check that element parameters can be set and retrieved
         */
        ImmersedBoundaryElement<2,2>* p_element = p_mesh->GetElement(0);

        p_element->SetMembraneSpringConstant(3.5);
        p_element->SetMembraneRestLength(1.7);

        double rest_length = p_mesh->GetElementIteratorBegin()->GetMembraneRestLength();
        double spring_constant = p_mesh->GetElementIteratorBegin()->GetMembraneSpringConstant();

        TS_ASSERT_EQUALS(rest_length, 1.7);
        TS_ASSERT_EQUALS(spring_constant, 3.5);
    }

    void xTestSuperellipseGenerator() throw(Exception)
    {
        SuperellipseGenerator* p_gen = new SuperellipseGenerator(100, 0.2, 0.2, 0.6, 0.2, 0.2);
        std::vector<c_vector<double, 2> > locations = p_gen->GetPointsAsVectors();

        std::vector<Node<2>*> nodes;

        for(unsigned location = 0 ; location < locations.size() ; location++)
        {
            nodes.push_back(new Node<2>(location, locations[location], true));
        }

        // Create a vector of immersed boundary elements and create an element with the nodes above
        std::vector<ImmersedBoundaryElement<2,2>*> ib_element;
        ib_element.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));

        // Create a mesh with the nodes and elements vectors
        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_element, 128, 128);

        ImmersedBoundaryElement<2,2>* p_elem = p_mesh->GetElement(0u);
        p_elem->SetMembraneRestLength(0.0001);
        p_elem->SetMembraneSpringConstant(10000.0);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

//        MAKE_PTR(ImmersedBoundaryElasticityForce<2>, p_elas_force);
//        p_main_modifier->AddImmersedBoundaryForce(p_elas_force);

        // Set simulation properties
        simulator.SetOutputDirectory("IB/TestSuperellipseGenerator");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(0.02);

        // Run the simulation
        simulator.Solve();
    }

    void xTestPalisadeGenerator() throw(Exception)
    {
        unsigned num_cells_wide    = 5;
        unsigned nodes_per_cell    = 100;
        double   ellipse_exponent  = 0.2;
        double   cell_aspect_ratio = 2.0;
        double   random_y_mult     = 0.15;
        bool     membrane          = true;

        ImmersedBoundaryPalisadeMeshGenerator gen(num_cells_wide,
                                                  nodes_per_cell,
                                                  ellipse_exponent,
                                                  cell_aspect_ratio,
                                                  random_y_mult,
                                                  membrane);

        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        p_mesh->GetMembraneElement()->SetMembraneSpringConstant(100000.0);
        p_mesh->GetMembraneElement()->SetMembraneRestLength(0.0001);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        // Set simulation properties
        simulator.SetOutputDirectory("IB/TestPalisadeGenerator");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(0.02);

        // Run the simulation
        //simulator.Solve();

    }

    void TestImmersedBoundarySimpleSimulation() throw(Exception)
    {
        // Load immersed boundary mesh
        ImmersedBoundaryMeshReader<2,2> mesh_reader("projects/ImmersedBoundary/test/mesh/ib_square_16");

        // Construct the immersed boundary mesh from the mesh reader
        ImmersedBoundaryMesh<2,2> mesh;

        mesh.ConstructFromMeshReader(mesh_reader);

        mesh.SetNumGridPtsX(64);
        mesh.SetNumGridPtsY(64);

        ImmersedBoundaryElement<2,2>* p_elem = mesh.GetElement(0u);
        p_elem->SetMembraneRestLength(0.005);
        p_elem->SetMembraneSpringConstant(1000.0);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

        PRINT_VARIABLE(cells.size());

        ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        MAKE_PTR(ImmersedBoundaryElasticityForce<2>, p_elas_force);
        p_main_modifier->AddImmersedBoundaryForce(p_elas_force);

        // Set simulation properties
        simulator.SetOutputDirectory("IB/TestImmersedBoundary");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(0.02);

        // Run the simulation
        simulator.Solve();
    }
};