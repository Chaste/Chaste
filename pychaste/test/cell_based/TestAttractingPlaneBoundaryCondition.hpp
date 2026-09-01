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

#ifndef TESTATTRACTINGPLANEBOUNDARYCONDITION_HPP_
#define TESTATTRACTINGPLANEBOUNDARYCONDITION_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "AttractingPlaneBoundaryCondition.hpp"

#include "ArchiveLocationInfo.hpp"
#include "ArchiveOpener.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileFinder.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Tests for AttractingPlaneBoundaryCondition, focusing on the
 * checkpoint (save/load) path that runs the constructor via load_construct_data.
 */
class TestAttractingPlaneBoundaryCondition : public AbstractCellBasedTestSuite
{
public:

    /**
     * Save and restore an AttractingPlaneBoundaryCondition through a Boost
     * archive. Deserialization reconstructs the object via load_construct_data,
     * which re-runs the constructor (now carrying the off-lattice check), so
     * this confirms the moved check does not throw on a valid population and
     * that every stored parameter survives the round-trip.
     */
    void TestArchivingOfAttractingPlaneBoundaryCondition()
    {
        EXIT_IF_PARALLEL; // Archiving of cell-based simulations is not supported in parallel.

        // Build a node-based (off-lattice) population.
        TrianglesMeshReader<2, 2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2, 2> generating_mesh;
        generating_mesh.ConstructFromMeshReader(mesh_reader);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> population(mesh, cells);

        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "AttractingPlaneBoundaryCondition.arch";
        ArchiveLocationInfo::SetMeshFilename("AttractingPlaneBoundaryCondition");

        c_vector<double, 2> point = zero_vector<double>(2);
        point[0] = 0.25;
        c_vector<double, 2> normal = unit_vector<double>(2, 1); // (0, 1)

        {
            // Constructing on an off-lattice population must not throw.
            AttractingPlaneBoundaryCondition<2> boundary_condition(&population, point, normal);
            boundary_condition.SetUseJiggledNodesOnPlane(true);
            boundary_condition.SetAttractionThreshold(0.5);

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Serialize via a base-class pointer, as a simulation would.
            AbstractCellPopulationBoundaryCondition<2>* const p_boundary_condition = &boundary_condition;
            (*p_arch) << p_boundary_condition;
        }

        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCellPopulationBoundaryCondition<2, 2>* p_boundary_condition;

            // Restore from the archive: this calls the constructor again via
            // load_construct_data and must not throw.
            (*p_arch) >> p_boundary_condition;

            auto p_restored = static_cast<AttractingPlaneBoundaryCondition<2>*>(p_boundary_condition);

            // The point and normal are reconstructed by load_construct_data.
            TS_ASSERT_DELTA(p_restored->rGetPointOnPlane()[0], 0.25, 1e-6);
            TS_ASSERT_DELTA(p_restored->rGetPointOnPlane()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_restored->rGetNormalToPlane()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_restored->rGetNormalToPlane()[1], 1.0, 1e-6);

            // The jiggle flag and attraction threshold are restored by serialize().
            TS_ASSERT_EQUALS(p_restored->GetUseJiggledNodesOnPlane(), true);
            TS_ASSERT_DELTA(p_restored->GetAttractionThreshold(), 0.5, 1e-6);

            // Tidy up. load_construct_data deserialized a fresh population onto
            // the heap, so delete it (via the public accessor) along with the
            // restored boundary condition.
            delete p_restored->GetCellPopulation();
            delete p_boundary_condition;
        }
    }

    /**
     * The off-lattice check now lives in the constructor, so constructing on a
     * lattice-based population should throw immediately rather than on the first
     * call to ImposeBoundaryCondition().
     */
    void TestConstructorRejectsOnLatticePopulation()
    {
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        boost::shared_ptr<PottsMesh<2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        PottsBasedCellPopulation<2> potts_population(*p_mesh, cells);

        c_vector<double, 2> point = zero_vector<double>(2);
        c_vector<double, 2> normal = unit_vector<double>(2, 0); // (1, 0)

        TS_ASSERT_THROWS_THIS(
            AttractingPlaneBoundaryCondition<2> boundary_condition(&potts_population, point, normal),
            "AttractingPlaneBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }
};

#endif // TESTATTRACTINGPLANEBOUNDARYCONDITION_HPP_
