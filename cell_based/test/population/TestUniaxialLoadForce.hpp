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

#ifndef TESTUNIAXIALLOADFORCE_HPP_
#define TESTUNIAXIALLOADFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include <fstream>
#include <string>
#include <vector>

#include "OutputFileHandler.hpp"
#include "UniaxialLoadForce.hpp"

#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "SemBasedCellPopulation.hpp"
#include "SemForce.hpp"
#include "SemMesh.hpp"
#include "SemSingleElementMeshGenerator.hpp"
#include "SemSphericalElementMeshGenerator.hpp"

#include "ForwardEulerNumericalMethod.hpp"
#include "OffLatticeSimulation.hpp"

#include "AbstractCellBasedTestSuite.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"


class TestUniaxialLoadForce : public AbstractCellBasedTestSuite
{
private:

    /**
     * @param rMesh the mesh whose nodes carry the applied forces
     * @return the sum of the applied force over every node
     */
    template <unsigned DIM>
    c_vector<double, DIM> SumOfAppliedForces(SemMesh<DIM>& rMesh)
    {
        c_vector<double, DIM> total = zero_vector<double>(DIM);

        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            total += rMesh.GetNode(i)->rGetAppliedForce();
        }

        return total;
    }

    /**
     * @param numCells the number of cells to generate
     * @return that many cells, each with no cell cycle model
     */
    template <unsigned DIM>
    std::vector<CellPtr> CreateCells(unsigned numCells)
    {
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, DIM> cells_generator;
        cells_generator.GenerateBasicRandom(cells, numCells);
        return cells;
    }

public:

    /**
     * A 3x3x3 grid of unit scale has nodes in three layers a third of a unit apart, so a slab half
     * a node spacing thick grips exactly the nine nodes of the outermost layer at each end.
     */
    void TestSlabIdentificationOnAKnownGrid()
    {
        SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 1.0);
        auto p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        UniaxialLoadForce<3> force;
        force.SetLoad(9.0);
        force.SetSlabThickness(0.5 / 3.0);

        // The slabs are not identified until the force is first used
        TS_ASSERT_EQUALS(force.GetUpperSlabNodeCount(), 0u);
        TS_ASSERT_EQUALS(force.GetLowerSlabNodeCount(), 0u);

        force.AddForceContribution(population);

        TS_ASSERT_EQUALS(force.GetUpperSlabNodeCount(), 9u);
        TS_ASSERT_EQUALS(force.GetLowerSlabNodeCount(), 9u);

        // The loading axis defaults to the last one, and the other settings round-trip
        TS_ASSERT_EQUALS(force.GetLoadingAxis(), 2u);
        TS_ASSERT_DELTA(force.GetLoad(), 9.0, 1e-12);
        TS_ASSERT_DELTA(force.GetSlabThickness(), 0.5 / 3.0, 1e-12);
        TS_ASSERT(force.GetIsCompressive());
    }

    /**
     * Each slab carries a total load of exactly the requested value, shared equally between its
     * nodes, and nodes gripped by neither slab are left alone.
     */
    void TestAppliedForcesAreExactAndSumToZero()
    {
        SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 1.0);
        auto p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        // Nine nodes per slab and a load of nine gives a per-node force of exactly one
        UniaxialLoadForce<3> force;
        force.SetLoad(9.0);
        force.SetSlabThickness(0.5 / 3.0);
        force.AddForceContribution(population);

        unsigned num_pushed_down = 0u;
        unsigned num_pushed_up = 0u;
        unsigned num_unloaded = 0u;

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            const c_vector<double, 3>& r_force = p_mesh->GetNode(i)->rGetAppliedForce();
            const double z = p_mesh->GetNode(i)->rGetLocation()[2];

            // Nothing is applied transverse to the loading axis
            TS_ASSERT_DELTA(r_force[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(r_force[1], 0.0, 1e-12);

            if (z > 0.5)
            {
                TS_ASSERT_DELTA(r_force[2], -1.0, 1e-12);
                num_pushed_down++;
            }
            else if (z < 1e-6)
            {
                TS_ASSERT_DELTA(r_force[2], 1.0, 1e-12);
                num_pushed_up++;
            }
            else
            {
                TS_ASSERT_DELTA(r_force[2], 0.0, 1e-12);
                num_unloaded++;
            }
        }

        TS_ASSERT_EQUALS(num_pushed_down, 9u);
        TS_ASSERT_EQUALS(num_pushed_up, 9u);
        TS_ASSERT_EQUALS(num_unloaded, 9u);

        const c_vector<double, 3> total = SumOfAppliedForces<3>(*p_mesh);
        TS_ASSERT_DELTA(norm_2(total), 0.0, 1e-12);
    }

    /**
     * The guarantee that matters. A ball carved from a lattice does not in general have the same
     * number of nodes at each end, so the two slabs carry different per-node forces. Dividing each
     * slab's load by its own node count is what keeps the totals equal and opposite; a common
     * per-node force would leave a residual of one node's worth and the specimen would drift.
     */
    void TestUnequalSlabsStillSumToZero()
    {
        SemSphericalElementMeshGenerator<3> generator(200, 0.5);
        auto p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        UniaxialLoadForce<3> force;
        force.SetLoad(1.0);
        force.SetSlabThickness(generator.GetNodeSpacing());
        force.AddForceContribution(population);

        // The two ends of this particular cell hold different numbers of nodes
        TS_ASSERT_DIFFERS(force.GetUpperSlabNodeCount(), force.GetLowerSlabNodeCount());
        TS_ASSERT_LESS_THAN(0u, force.GetUpperSlabNodeCount());
        TS_ASSERT_LESS_THAN(0u, force.GetLowerSlabNodeCount());

        // Each slab nonetheless carries a total of exactly the requested load
        double upper_total = 0.0;
        double lower_total = 0.0;
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            const double applied = p_mesh->GetNode(i)->rGetAppliedForce()[2];
            if (applied < 0.0)
            {
                upper_total += applied;
            }
            else
            {
                lower_total += applied;
            }
        }

        TS_ASSERT_DELTA(upper_total, -1.0, 1e-12);
        TS_ASSERT_DELTA(lower_total, 1.0, 1e-12);

        const c_vector<double, 3> total = SumOfAppliedForces<3>(*p_mesh);
        TS_ASSERT_DELTA(norm_2(total), 0.0, 1e-12);
    }

    /**
     * Tension is compression with both loads reversed, so it must produce exactly the opposite
     * force on every node.
     */
    void TestTensionActsOppositeToCompression()
    {
        SemSingleElementMeshGenerator<3> compressive_generator({ 3, 3, 3 }, 1.0);
        auto p_compressive_mesh = compressive_generator.GetMesh();
        std::vector<CellPtr> compressive_cells = CreateCells<3>(p_compressive_mesh->GetNumElements());
        SemBasedCellPopulation<3> compressive_population(*p_compressive_mesh, compressive_cells);

        UniaxialLoadForce<3> compressive_force;
        compressive_force.SetLoad(9.0);
        compressive_force.SetSlabThickness(0.5 / 3.0);
        compressive_force.AddForceContribution(compressive_population);

        SemSingleElementMeshGenerator<3> tensile_generator({ 3, 3, 3 }, 1.0);
        auto p_tensile_mesh = tensile_generator.GetMesh();
        std::vector<CellPtr> tensile_cells = CreateCells<3>(p_tensile_mesh->GetNumElements());
        SemBasedCellPopulation<3> tensile_population(*p_tensile_mesh, tensile_cells);

        UniaxialLoadForce<3> tensile_force;
        tensile_force.SetLoad(9.0);
        tensile_force.SetSlabThickness(0.5 / 3.0);
        tensile_force.SetIsCompressive(false);
        TS_ASSERT(!tensile_force.GetIsCompressive());
        tensile_force.AddForceContribution(tensile_population);

        for (unsigned i = 0; i < p_compressive_mesh->GetNumNodes(); ++i)
        {
            TS_ASSERT_DELTA(p_tensile_mesh->GetNode(i)->rGetAppliedForce()[2],
                            -p_compressive_mesh->GetNode(i)->rGetAppliedForce()[2], 1e-12);
        }
    }

    /**
     * The slabs are gripped once and then held, as a rheometer plate grips a specimen, so moving
     * the nodes afterwards must not hand the load to a different set of them.
     */
    void TestSlabMembershipIsFixedAfterNodesMove()
    {
        SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 1.0);
        auto p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        UniaxialLoadForce<3> force;
        force.SetLoad(9.0);
        force.SetSlabThickness(0.5 / 3.0);
        force.AddForceContribution(population);

        // Record which nodes were gripped, then collapse the specimen onto a plane so that every
        // node would qualify for both slabs were membership recomputed
        std::vector<double> original_forces;
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            original_forces.push_back(p_mesh->GetNode(i)->rGetAppliedForce()[2]);
            p_mesh->GetNode(i)->rGetModifiableLocation()[2] = 0.0;
            p_mesh->GetNode(i)->ClearAppliedForce();
        }

        force.AddForceContribution(population);

        TS_ASSERT_EQUALS(force.GetUpperSlabNodeCount(), 9u);
        TS_ASSERT_EQUALS(force.GetLowerSlabNodeCount(), 9u);

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetAppliedForce()[2], original_forces[i], 1e-12);
        }
    }

    /**
     * End to end: a cell held together by SemForce and compressed by this force gets shorter along
     * the loading axis, and does not wander off while doing so.
     */
    void TestCompressionShortensTheSpecimen()
    {
        // Mesh, box collection and force parameters follow the proven 3D case of
        // TestRunningSemBasedSimulationsTutorial, so that only the load is under test here
        SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 0.5);
        auto p_mesh = generator.GetMesh();
        const double node_spacing = 0.5 / 3.0;

        const double interaction_cutoff = 0.25;
        c_vector<double, 6> box_domain;
        box_domain[0] = -1.0;
        box_domain[1] = 2.0;
        box_domain[2] = -1.0;
        box_domain[3] = 2.0;
        box_domain[4] = -1.0;
        box_domain[5] = 2.0;
        p_mesh->SetUpBoxCollection(interaction_cutoff, box_domain);

        const unsigned num_nodes = p_mesh->GetNumNodes();
        const double cell_radius = 0.25;
        const double kappa0 = 20.0;
        const double rho = 5.0;
        const double packing = 1.0;
        const double eta0 = 1.0 / static_cast<double>(num_nodes);

        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        const double initial_height = p_mesh->GetWidth(2u);
        const c_vector<double, 3> initial_centroid = p_mesh->GetCentroidOfElement(0u);

        OffLatticeSimulation<3> simulator(population);
        simulator.SetOutputDirectory("TestUniaxialLoadForceCompression");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(1.0);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<3> >());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

        MAKE_PTR(SemForce<3>, p_sem_force);
        p_sem_force->SetIntraScalingFactor(rho);
        const SemNScaledParameters nscaled
            = p_sem_force->ApplyNScaledIntraParameters(num_nodes, cell_radius, kappa0, 0.0, packing, eta0);
        p_sem_force->SetIntraCutOffDistance(interaction_cutoff);
        population.SetDampingConstantNormal(nscaled.DampingConstant);
        simulator.AddForce(p_sem_force);

        MAKE_PTR(UniaxialLoadForce<3>, p_load);
        p_load->SetLoad(0.5);
        p_load->SetSlabThickness(0.5 * node_spacing);
        simulator.AddForce(p_load);

        simulator.Solve();

        // The specimen is shorter, by a strain of a few percent: enough to be unambiguous, small
        // enough that the response is still that of the intact cell rather than of a collapsed one
        const double final_height = p_mesh->GetWidth(2u);
        const double strain = (initial_height - final_height) / initial_height;
        TS_ASSERT_LESS_THAN(final_height, initial_height);
        TS_ASSERT_DELTA(strain, 0.033, 0.015);

        // Equal and opposite loads leave the centre of mass where it was
        const c_vector<double, 3> final_centroid = p_mesh->GetCentroidOfElement(0u);
        TS_ASSERT_DELTA(norm_2(final_centroid - initial_centroid), 0.0, 1e-6);
    }

    /**
     * Nothing about this force is specific to subcellular element populations. Applied to a
     * node-based population, where one node is one cell, the same protocol compresses a tissue
     * rather than a cell.
     */
    void TestWorksWithANodeBasedPopulation()
    {
        // Two layers of four cells, a unit apart
        std::vector<Node<3>*> nodes;
        unsigned index = 0u;
        for (unsigned k = 0; k < 2; ++k)
        {
            for (unsigned j = 0; j < 2; ++j)
            {
                for (unsigned i = 0; i < 2; ++i)
                {
                    nodes.push_back(new Node<3>(index++, false,
                                                static_cast<double>(i),
                                                static_cast<double>(j),
                                                static_cast<double>(k)));
                }
            }
        }

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells = CreateCells<3>(mesh.GetNumNodes());
        NodeBasedCellPopulation<3> population(mesh, cells);

        UniaxialLoadForce<3> force;
        force.SetLoad(4.0);
        force.SetSlabThickness(0.4);
        force.AddForceContribution(population);

        TS_ASSERT_EQUALS(force.GetUpperSlabNodeCount(), 4u);
        TS_ASSERT_EQUALS(force.GetLowerSlabNodeCount(), 4u);

        c_vector<double, 3> total = zero_vector<double>(3);
        for (unsigned i = 0; i < mesh.GetNumNodes(); ++i)
        {
            const c_vector<double, 3>& r_force = population.GetNode(i)->rGetAppliedForce();
            total += r_force;

            // A load of four over four nodes per slab is one apiece
            TS_ASSERT_DELTA(fabs(r_force[2]), 1.0, 1e-12);
        }

        TS_ASSERT_DELTA(norm_2(total), 0.0, 1e-12);

        for (unsigned i = 0; i < nodes.size(); ++i)
        {
            delete nodes[i];
        }
    }

    /**
     * The slabs are identified once and then held fixed, so a restart must recover not just the
     * loading parameters but the slab membership itself: re-identifying it after the specimen has
     * deformed would pick a different set of nodes.
     */
    void TestArchiving()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "UniaxialLoadForce.arch";

        unsigned upper_slab_count = 0u;
        unsigned lower_slab_count = 0u;

        {
            SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 1.0);
            auto p_mesh = generator.GetMesh();

            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 3> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            SemBasedCellPopulation<3> population(*p_mesh, cells);

            UniaxialLoadForce<3> force;
            force.SetLoad(2.5);
            force.SetLoadingAxis(1u);
            force.SetSlabThickness(0.2);
            force.SetIsCompressive(false);

            // Applying the force identifies the slabs, so they become part of the archived state
            force.AddForceContribution(population);
            upper_slab_count = force.GetUpperSlabNodeCount();
            lower_slab_count = force.GetLowerSlabNodeCount();
            TS_ASSERT_LESS_THAN(0u, upper_slab_count);
            TS_ASSERT_LESS_THAN(0u, lower_slab_count);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer to most abstract class possible
            AbstractForce<3>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<3>* p_force;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_force;

            UniaxialLoadForce<3>* p_loaded = static_cast<UniaxialLoadForce<3>*>(p_force);
            TS_ASSERT_DELTA(p_loaded->GetLoad(), 2.5, 1e-12);
            TS_ASSERT_EQUALS(p_loaded->GetLoadingAxis(), 1u);
            TS_ASSERT_DELTA(p_loaded->GetSlabThickness(), 0.2, 1e-12);
            TS_ASSERT_EQUALS(p_loaded->GetIsCompressive(), false);

            // The slab membership survives, so the restarted force loads the same nodes
            TS_ASSERT_EQUALS(p_loaded->GetUpperSlabNodeCount(), upper_slab_count);
            TS_ASSERT_EQUALS(p_loaded->GetLowerSlabNodeCount(), lower_slab_count);

            delete p_force;
        }
    }

    void TestExceptions()
    {
        UniaxialLoadForce<3> force;

        TS_ASSERT_THROWS_THIS(force.SetLoad(0.0),
            "UniaxialLoadForce: load must be positive");
        TS_ASSERT_THROWS_THIS(force.SetLoad(-1.0),
            "UniaxialLoadForce: load must be positive");

        TS_ASSERT_THROWS_THIS(force.SetLoadingAxis(3u),
            "UniaxialLoadForce: loadingAxis must be less than DIM");
        TS_ASSERT_THROWS_NOTHING(force.SetLoadingAxis(2u));

        TS_ASSERT_THROWS_THIS(force.SetSlabThickness(0.0),
            "UniaxialLoadForce: slabThickness must be positive");
        TS_ASSERT_THROWS_THIS(force.SetSlabThickness(-1.0),
            "UniaxialLoadForce: slabThickness must be positive");

        SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 1.0);
        auto p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        // Neither the load nor the slab thickness has been set yet
        UniaxialLoadForce<3> unconfigured_force;
        TS_ASSERT_THROWS_THIS(unconfigured_force.AddForceContribution(population),
            "UniaxialLoadForce: SetLoad() must be called before the force is used");

        unconfigured_force.SetLoad(1.0);
        TS_ASSERT_THROWS_THIS(unconfigured_force.AddForceContribution(population),
            "UniaxialLoadForce: SetSlabThickness() must be called before the force is used");

        // The specimen is two thirds of a unit deep, so slabs half a unit thick would overlap
        UniaxialLoadForce<3> overlapping_force;
        overlapping_force.SetLoad(1.0);
        overlapping_force.SetSlabThickness(0.5);
        TS_ASSERT_THROWS_THIS(overlapping_force.AddForceContribution(population),
            "UniaxialLoadForce: slabThickness is more than half the extent of the population "
            "along the loading axis, so the two slabs would overlap");
    }
};

#endif /*TESTUNIAXIALLOADFORCE_HPP_*/
