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

// Needed for the test environment
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "OffLatticeSimulation.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"
#include "OutputFileHandler.hpp"

#include "Debug.hpp"
#include "Timer.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestNumericsPaperSimulations : public AbstractCellBasedTestSuite
{
public:

    void TestEllipseRelaxing() throw(Exception)
    {
        /*
         * 1: num nodes
         * 2: superellipse exponent
         * 3: cell width
         * 4: cell height
         * 5: bottom left x
         * 6: bottom left y
         */
        SuperellipseGenerator* p_gen = new SuperellipseGenerator(128, 1.0, 0.4, 0.6, 0.3, 0.2);
        std::vector<c_vector<double, 2> > locations = p_gen->GetPointsAsVectors();

        std::vector<Node<2>* > nodes;
        std::vector<ImmersedBoundaryElement<2,2>* > elements;

        for(unsigned location = 0 ; location < locations.size() ; location++)
        {
            nodes.push_back(new Node<2>(location, locations[location], true));
        }

        elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));
        elements[0]->rGetCornerNodes().push_back(nodes[0]);
        elements[0]->rGetCornerNodes().push_back(nodes[1]);
        elements[0]->rGetCornerNodes().push_back(nodes[2]);
        elements[0]->rGetCornerNodes().push_back(nodes[3]);

        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2,2>(nodes, elements);
        p_mesh->SetNumGridPtsXAndY(32);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(false);

        OffLatticeSimulation<2> simulator(cell_population);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        MAKE_PTR_ARGS(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force, (cell_population));
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetSpringConstant(1e8);

        std::string output_directory = "numerics_paper/ellipse_relaxing";
        simulator.SetOutputDirectory(output_directory);

        // Write the state of the immersed boundary mesh to file
        ImmersedBoundaryMeshWriter<2,2> mesh_at_start(output_directory, "example_simulation_mesh_at_start");
        ImmersedBoundaryMeshWriter<2,2> mesh_at_end(output_directory, "example_simulation_mesh_at_end");

        OutputFileHandler results_handler(output_directory, false);
        out_stream results_file = results_handler.OpenOutputFile("example_simulation_esf.dat");

        // Output summary statistics to results file
        (*results_file) << "time,esf\n";
        (*results_file) << 0.0 << "," << p_mesh->GetElongationShapeFactorOfElement(0) << "\n";

        // Set simulation properties
        double dt = 0.05;
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(1);

        for (unsigned i=0 ; i< 100 ; i++)
        {
            double new_end_time = dt * (1.0 + i);

            simulator.SetEndTime(new_end_time);
            simulator.Solve();

            (*results_file) << new_end_time << "," << p_mesh->GetElongationShapeFactorOfElement(0) << "\n";

            if(i==0)
            {
                mesh_at_start.WriteFilesUsingMesh(*p_mesh);
            }
        }

        mesh_at_end.WriteFilesUsingMesh(*p_mesh);

        // Tidy up
        results_file->close();
        delete(p_mesh);
    }
};