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

#ifndef TESTPARALLELCELLSGENERATOR_HPP_
#define TESTPARALLELCELLSGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "ParallelCellsGenerator.hpp"
#include "Cell.hpp"
#include "NodesOnlyMesh.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "CellPropertyRegistry.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"

class TestParallelCellsGenerator : public AbstractCellBasedTestSuite
{
public:

    void TestConstruct1dParallelPopulation()
    {
        NodesOnlyMesh<1> mesh;
        mesh.SetMaximumInteractionDistance(0.5);

        std::vector<CellPtr> cells;

        ParallelCellsGenerator<FixedG1GenerationalCellCycleModel, 1> generator;
        c_vector<double, 2> bounding_box = generator.GetArchiveBoundingBox("cell_based/test/data/TestParallelConstruction/Population1d.dat");

        TS_ASSERT_DELTA(bounding_box[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(bounding_box[1], 1.0, 1e-4);

        // Check an assertion which tests whether the space dimensions of the mesh and archive file are consistent.
        TS_ASSERT_THROWS_THIS(generator.GetArchiveBoundingBox("cell_based/test/data/TestParallelConstruction/Population2d.dat"), "Space dimension of ParallelCellsGenerator and archive file do not match");

        generator.GenerateParallelCells("cell_based/test/data/TestParallelConstruction/Population1d.dat", cells, mesh, CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

        unsigned num_nodes = mesh.GetNumNodes();
        unsigned total_nodes;

        MPI_Allreduce(&num_nodes, &total_nodes, 1, MPI_UNSIGNED, MPI_SUM, PetscTools::GetWorld());
        TS_ASSERT_EQUALS(total_nodes, 2u);

        unsigned num_cells = cells.size();
        unsigned total_cells;

        MPI_Allreduce(&num_cells, &total_cells, 1, MPI_UNSIGNED, MPI_SUM, PetscTools::GetWorld());
        TS_ASSERT_EQUALS(total_cells, 2u);

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT(cells[i]->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>());
        }
    }

    void TestConstruct2dParallelPopulation()
    {
        NodesOnlyMesh<2> mesh;
        mesh.SetMaximumInteractionDistance(0.5);

        std::vector<CellPtr> cells;

        ParallelCellsGenerator<UniformG1GenerationalCellCycleModel, 2> generator;
        c_vector<double, 4> bounding_box = generator.GetArchiveBoundingBox("cell_based/test/data/TestParallelConstruction/Population2d.dat");

        TS_ASSERT_DELTA(bounding_box[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(bounding_box[1], 1.0, 1e-4);
        TS_ASSERT_DELTA(bounding_box[2], 0.0, 1e-4);
        TS_ASSERT_DELTA(bounding_box[3], 1.0, 1e-4);

        generator.GenerateParallelCells("cell_based/test/data/TestParallelConstruction/Population2d.dat", cells, mesh, CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        unsigned num_nodes = mesh.GetNumNodes();
        unsigned total_nodes;

        MPI_Allreduce(&num_nodes, &total_nodes, 1, MPI_UNSIGNED, MPI_SUM, PetscTools::GetWorld());
        TS_ASSERT_EQUALS(total_nodes, 4u);

        unsigned num_cells = cells.size();
        unsigned total_cells;

        MPI_Allreduce(&num_cells, &total_cells, 1, MPI_UNSIGNED, MPI_SUM, PetscTools::GetWorld());
        TS_ASSERT_EQUALS(total_cells, 4u);

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT(cells[i]->GetCellProliferativeType()->IsType<TransitCellProliferativeType>());
        }
    }

    void TestConstruct3dParallelPopulation()
    {
        NodesOnlyMesh<3> mesh;
        mesh.SetMaximumInteractionDistance(0.5);

        std::vector<CellPtr> cells;

        ParallelCellsGenerator<UniformG1GenerationalCellCycleModel, 3> generator;
        c_vector<double, 6> bounding_box = generator.GetArchiveBoundingBox("cell_based/test/data/TestParallelConstruction/Population3d.dat");

        TS_ASSERT_DELTA(bounding_box[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(bounding_box[1], 1.0, 1e-4);
        TS_ASSERT_DELTA(bounding_box[2], 0.0, 1e-4);
        TS_ASSERT_DELTA(bounding_box[3], 1.0, 1e-4);
        TS_ASSERT_DELTA(bounding_box[4], 0.0, 1e-4);
        TS_ASSERT_DELTA(bounding_box[5], 1.0, 1e-4);

        generator.GenerateParallelCells("cell_based/test/data/TestParallelConstruction/Population3d.dat", cells, mesh, CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        unsigned num_nodes = mesh.GetNumNodes();
        unsigned total_nodes;

        MPI_Allreduce(&num_nodes, &total_nodes, 1, MPI_UNSIGNED, MPI_SUM, PetscTools::GetWorld());
        TS_ASSERT_EQUALS(total_nodes, 8u);

        unsigned num_cells = cells.size();
        unsigned total_cells;

        MPI_Allreduce(&num_cells, &total_cells, 1, MPI_UNSIGNED, MPI_SUM, PetscTools::GetWorld());
        TS_ASSERT_EQUALS(total_cells, 8u);

        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT(cells[i]->GetCellProliferativeType()->IsType<TransitCellProliferativeType>());
        }
    }
};

#endif /*TESTPARALLELCELLSGENERATOR_HPP_*/
