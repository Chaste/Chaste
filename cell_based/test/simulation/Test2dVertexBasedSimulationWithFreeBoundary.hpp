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

#ifndef TEST2DVERTEXBASEDSIMULATIONWITHFREEBOUNDARY_HPP_
#define TEST2DVERTEXBASEDSIMULATIONWITHFREEBOUNDARY_HPP_

#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "FarhadifarForce.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"
#include "CellsGenerator.hpp"
#include "VoronoiVertexMeshGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"

/**
 * This class consists of a single test, in which a 2D model
 * of a growing monolayer of cells is simulated for a fixed
 * period of time.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class Test2DVertexSimulationWithFreeBoundary : public AbstractCellBasedTestSuite
{
public:

    void Test2DFreeBoundaryVertexSimulationForProfiling()
    {
        // make the simulation
        std::vector<CellPtr> cells;

        MutableVertexMesh<2,2>* p_mesh;
        VoronoiVertexMeshGenerator mesh_generator = VoronoiVertexMeshGenerator(35,35,1,1.0);
        p_mesh = mesh_generator.GetMesh();
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;

        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        double initial_target_area = 1.0;

        for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            // target areas
            cell_iter->GetCellData()->SetItem("target area", initial_target_area);
        }

        OffLatticeSimulation<2> simulator(cell_population);

        // Make the Farhadifar force
        MAKE_PTR(FarhadifarForce<2>, p_force);

        // before passing the force to the simulation
        simulator.AddForce(p_force);

        // We need a FarhadifarType target area modifier
        MAKE_PTR(TargetAreaLinearGrowthModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.SetEndTime(20);
        simulator.SetDt(0.01);

        //    mpSimulator->SetSamplingTimestepMultiple( 100 );

        simulator.SetOutputDirectory("VertexSimulationWithFreeBoundary");

        simulator.Solve();
    }
};

#endif /*TEST2DVERTEXBASEDSIMULATIONWITHFREEBOUNDARY_HPP_*/
