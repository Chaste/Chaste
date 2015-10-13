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

// External library - not part of Chaste
#include <fftw3.h>

#include "OffLatticeSimulation.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "SuperellipseGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"

#include "Debug.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestShortSingleCellSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestShortSingleCellSim() throw(Exception)
    {
        /*
         * Create an Immersed Boundary Mesh using a SuperellipseGenerator
         *
         * 1: Num nodes
         * 2: Superellipse exponent
         * 3: Width
         * 4: Height
         * 5: Bottom-left x
         * 6: Botton-left y
         */
        SuperellipseGenerator *p_gen = new SuperellipseGenerator(128, 0.2, 0.3, 0.6, 0.35, 0.2);

        // Generate a mesh using this superellipse
        std::vector<c_vector<double, 2> > locations = p_gen->GetPointsAsVectors();
        delete p_gen;

        std::vector<Node<2>* > nodes;
        for (unsigned node_idx = 0 ; node_idx < locations.size() ; node_idx++)
        {
            nodes.push_back(new Node<2>(node_idx, locations[node_idx], true));
        }

        std::vector<ImmersedBoundaryElement<2,2>* > elements;
        elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));

        ImmersedBoundaryMesh<2,2> mesh(nodes, elements, 2048, 2048);

        mesh.GetElement(0)->SetMembraneSpringConstant(1e7);
        mesh.GetElement(0)->SetMembraneRestLength(0.25 * mesh.GetCharacteristicNodeSpacing());

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        MAKE_PTR_ARGS(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force, (cell_population));
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);


        // Set simulation properties
        simulator.SetOutputDirectory("TestShortSingleCellSimulation");
        simulator.SetDt(0.05);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(0.1);
        simulator.Solve();
    }
};