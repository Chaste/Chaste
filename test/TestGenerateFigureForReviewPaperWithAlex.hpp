/*

Copyright (c) 2005-2016, University of Oxford.
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

// Needed for the test environment
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// Includes from trunk
#include "CellId.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"

// Includes from Immersed Boundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestGenerateFigureForReviewPaperWithAlex : public AbstractCellBasedTestSuite
{
public:

    /*
     * Function to return a vector of node locations specifying a hexagon with origin at (0, 0) with vertices on the
     * unit circle.
     */
    std::vector<c_vector<double, 2> > GetUnitHexagon(unsigned numPtsPerSide)
    {
        std::vector<c_vector<double, 2> > locations(numPtsPerSide * 6);

        // Find locations of the six vertices (and the seventh is the same as the first)
        std::vector<c_vector<double, 2> > vertices(7);
        double sixty_degrees = M_PI / 3.0;
        for (unsigned vertex = 0; vertex < vertices.size(); vertex++)
        {
            vertices[vertex][0] = cos((double)vertex * sixty_degrees);
            vertices[vertex][1] = sin((double)vertex * sixty_degrees);
        }

        // Between each pair of vertices, fill in the specified number of locations evenly spaced along the edge
        for (unsigned vertex = 0; vertex < 6; vertex++)
        {
            c_vector<double, 2> this_vertex = vertices[vertex];
            c_vector<double, 2> next_vertex = vertices[vertex + 1];
            c_vector<double, 2> vec_between = next_vertex - this_vertex;

            for (unsigned i = 0; i < numPtsPerSide; i++)
            {
                locations[vertex * numPtsPerSide + i] = this_vertex + i * vec_between / (double)numPtsPerSide;
            }
        }
        return locations;
    }

    void TestHexagonalPacking() throw(Exception)
    {
        /*
         * Create a hexagonal packing of elements of a fixed radius
         */

        // Specify the number of elements in each direction and their radius
        unsigned num_x = 5;
        unsigned num_y = 4;
        double rad = 0.11;

        // We need to calculate the required centre of each hexagon
        std::vector<c_vector<double, 2> > offsets(num_x * num_y);

        // Specify the centre of the bottom-left most hexagon
        c_vector<double, 2> global_offset;
        global_offset[0] = 0.15;
        global_offset[1] = 0.15;

        // Calculate the centres
        for (unsigned x = 0; x < num_x; x++)
        {
            for (unsigned y = 0; y < num_y; y++)
            {
                unsigned idx = x * num_y + y;

                offsets[idx][0] = global_offset[0] + 1.5 * rad * (double)x;
                offsets[idx][1] = global_offset[1] + 0.5 * sqrt(3) * rad * (double)(2 * y + (x % 2 == 1));
            }
        }

        // Get locations for the reference hexagon
        std::vector<c_vector<double, 2> > node_locations = GetUnitHexagon(10);

        std::vector<Node<2>*> nodes;
        std::vector<ImmersedBoundaryElement<2,2>*> elems;

        // For each calculated centre, create the nodes representing each location around that hexagon
        for (unsigned offset = 0; offset < offsets.size(); offset++)
        {
            std::vector<Node<2>*> nodes_this_elem;

            for (unsigned location = 0; location < node_locations.size(); location++)
            {
                unsigned index = offset * node_locations.size() + location;
                Node<2>* p_node = new Node<2>(index, offsets[offset] + 0.95 * rad * node_locations[location], true);

                nodes_this_elem.push_back(p_node);
                nodes.push_back(p_node);
            }

            ImmersedBoundaryElement<2,2>* p_elem = new ImmersedBoundaryElement<2,2>(offset, nodes_this_elem);
            elems.push_back(p_elem);
        }

        // Create the mesh, cells, cell population, and simulation
        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2,2>(nodes, elems);
        p_mesh->SetNumGridPtsXAndY(64);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(true);

        OffLatticeSimulation<2> simulator(cell_population);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        MAKE_PTR(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force);
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetSpringConstant(0.5 * 1e7);

        // Set output location
        std::string output_directory = "TestGenerateFigureForReviewPaperWithAlex";
        simulator.SetOutputDirectory(output_directory);

        // Prepare mesh writers to export the .node and .grid files needed to draw quiver plots
        ImmersedBoundaryMeshWriter<2,2> mesh_at_1(output_directory, "mesh_at_1");
        ImmersedBoundaryMeshWriter<2,2> mesh_at_2(output_directory, "mesh_at_2");

        // Set simulation properties
        simulator.SetOutputDirectory(output_directory);
        double dt = 0.05;
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(100);

        /*
         * Run a first simulation to allow cells to relax to an equilibrium state, then write the mesh files
         */
        simulator.SetEndTime(2500.0 * dt);
        simulator.Solve();

        mesh_at_1.WriteFilesUsingMesh(*p_mesh);

        /*
         * Add a fluid sink to one of the middle cells to make it reduce in size.  Run the simulation for a further
         * duration and write the new mesh files.
         */
        std::vector<FluidSource<2>*>& sources = p_mesh->rGetElementFluidSources();
        sources[10]->SetStrength(-2.0 * 1e4);

        simulator.SetEndTime(7000.0 * dt);
        simulator.Solve();

        mesh_at_2.WriteFilesUsingMesh(*p_mesh);

        // Tidy up by deleting the mesh.  This will free memory used for nodes and elements.
        delete(p_mesh);
    }
};
