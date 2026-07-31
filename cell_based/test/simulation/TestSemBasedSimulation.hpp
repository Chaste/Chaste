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

#ifndef TESTSEMBASEDSIMULATION_HPP_
#define TESTSEMBASEDSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <cmath>

#include "CheckpointArchiveTypes.hpp"
#include "OffLatticeSimulation.hpp"
#include "SemSingleElementMeshGenerator.hpp"
#include "SemMultiElementMeshGenerator.hpp"
#include "SemSphericalElementMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "ToroidalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "ChemotacticForce.hpp"
#include "RandomCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulationWithMyStoppingEvent.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "CellIdWriter.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "SemBasedCellPopulation.hpp"
#include "SemForce.hpp"
#include "SemGaussianRandomForce.hpp"
#include "SemSpatiallyCorrelatedRandomForce.hpp"
#include "SemRegionalForce.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeLocationWriter.hpp"
#include "NodeRegionPointDataWriter.hpp"
#include "ElementIdNodePointDataWriter.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestSemBasedSimulation : public AbstractCellBasedWithTimingsTestSuite
{
private:

    /**
     * The RMS distance of the population's nodes from their centroid. Used in preference to a
     * bounding box because a bounding box is fixed by whichever two nodes happen to lie furthest
     * apart, and a close packing has a systematically different reach along different axes.
     *
     * @param rCellPopulation the population to measure
     *
     * @return the radius of gyration
     */
    double CalculateRadiusOfGyration(AbstractCellPopulation<3>& rCellPopulation)
    {
        c_vector<double, 3> centroid = zero_vector<double>(3);
        unsigned num_nodes = 0u;

        for (auto node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            centroid += node_iter->rGetLocation();
            num_nodes++;
        }
        assert(num_nodes > 0u);
        centroid /= static_cast<double>(num_nodes);

        double sum_of_squares = 0.0;
        for (auto node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            sum_of_squares += norm_2(node_iter->rGetLocation() - centroid)
                              * norm_2(node_iter->rGetLocation() - centroid);
        }

        return std::sqrt(sum_of_squares / static_cast<double>(num_nodes));
    }

public:

    void TestSemBasedSimulationExample2D()
    {
        SemSingleElementMeshGenerator<2> generator({3,3}, 0.5);
        auto p_mesh = generator.GetMesh();

        c_vector<double, 4> boxCollectionDomain{};
        boxCollectionDomain[0] = -1.0;
        boxCollectionDomain[1] =  2.0;
        boxCollectionDomain[2] = -1.0;
        boxCollectionDomain[3] =  2.0;

        p_mesh->SetUpBoxCollection(0.1, boxCollectionDomain);

        // Assertions
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 9);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.0);

        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 9);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1);


        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1);
        simulator.SetOutputDirectory("TestSemBasedSimulation");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(0.03);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

        // Create some force laws and pass them to the simulation
        MAKE_PTR(SemRegionalForce<2>, p_sem_force);
        simulator.AddForce(p_sem_force);

        MAKE_PTR(SemGaussianRandomForce<2>, p_random_force);
        p_random_force->SetDiffusionConstant(1e-6);
        simulator.AddForce(p_random_force);

        // Run the simulation
        simulator.Solve();
    }

    void TestSemBasedSimulationExample3D()
    {
        SemMultiElementMeshGenerator<3> generator({ 3, 3, 3 }, {1, 1, 1}, 0.5);
        auto p_mesh = generator.GetMesh();

        c_vector<double, 6> boxCollectionDomain{};
        boxCollectionDomain[0] = -1.0;
        boxCollectionDomain[1] =  2.0;
        boxCollectionDomain[2] = -1.0;
        boxCollectionDomain[3] =  2.0;
        boxCollectionDomain[4] = -1.0;
        boxCollectionDomain[5] =  2.0;

        const double interaction_cutoff = 0.25;
        p_mesh->SetUpBoxCollection(interaction_cutoff, boxCollectionDomain);

        // Assertions
        const unsigned expected_num_elements = 1u;
        const unsigned expected_num_nodes = expected_num_elements * 3u * 3u * 3u;
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), expected_num_elements);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), expected_num_nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddNodePointDataWriter<NodeRegionPointDataWriter>();
        cell_population.AddNodePointDataWriter<ElementIdNodePointDataWriter>();

        TS_ASSERT_EQUALS(cell_population.GetNumElements(), expected_num_elements);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), expected_num_nodes);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), expected_num_elements);


        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSemBasedSimulation3D");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(0.03);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<3>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

        // Create some force laws and pass them to the simulation.
        // R_cell = scaleFactor/2, packing=1 (regular cubic grid), kappa0 chosen to match well_depth~0.001.
        // The intra scaling also yields the consistent damping constant eta = eta0*N, which we
        // apply to the population so the Langevin dynamics are scaled as in Sandersius 2008 Section 2.
        MAKE_PTR(SemForce<3>, p_sem_force);
        const SemNScaledParameters intra_params = p_sem_force->ApplyNScaledIntraParameters(p_mesh->GetNumNodes(), 0.25, 20.0, 0.0, 1.0);
        p_sem_force->SetIntraCutOffDistance(interaction_cutoff);
        p_sem_force->ApplyNScaledInterParameters(p_mesh->GetNumNodes(), 0.25, 20.0, 0.0, 1.0);
        p_sem_force->SetInterCutOffDistance(interaction_cutoff);
        cell_population.SetDampingConstantNormal(intra_params.DampingConstant);
        simulator.AddForce(p_sem_force);

        MAKE_PTR(SemSpatiallyCorrelatedRandomForce<3>, p_random_force);
        p_random_force->SetDiffusionConstant(1 * 1e-5);
        p_random_force->SetCorrelationLength(p_sem_force->GetIntraEquilibriumDistance());
        p_random_force->SetLowerCorner({{-1.0, -1.0, -1.0}});
        p_random_force->SetUpperCorner({{3.0, 3.0, 3.0}});
        simulator.AddForce(p_random_force);

        // Run the simulation
        simulator.Solve();
    }

    /**
     * The spherical generator is the paper-faithful initial condition, but it is otherwise only
     * tested at mesh level. This case runs it through a population and a short relaxation, and in
     * doing so checks the generator's own node spacing against the paper's r_eq relation on a mesh
     * that actually has that spacing -- elsewhere that consistency is asserted only in the abstract
     * by TestSemParameterScaling.
     */
    void TestSemBasedSimulationWithSphericalElement()
    {
        const unsigned num_nodes = 200u;
        const double cell_radius = 0.25;

        // p = pi/(3*sqrt(2)), the FCC/HCP close packing density. The generator carves the cell out
        // of a perfect close packing, so this is the packing that applies -- not the p = 0.64 of a
        // random close packing that the paper's own equilibrated aggregates reach.
        const double packing_density = M_PI / (3.0 * std::sqrt(2.0));

        SemSphericalElementMeshGenerator<3> generator(num_nodes, cell_radius);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), num_nodes);

        /*
         * The generator's spacing comes from geometry (carve N sites out of a close packing, scale
         * so the outermost sits at cell_radius); the paper's r_eq comes from a packing argument.
         * They are independent derivations, and agree to about 2%.
         */
        const double paper_r_eq = 2.0 * cell_radius
            * std::pow(packing_density / static_cast<double>(num_nodes), 1.0 / 3.0);
        TS_ASSERT_DELTA(generator.GetNodeSpacing(), paper_r_eq, 0.02 * paper_r_eq);

        const double interaction_cutoff = 2.0 * generator.GetNodeSpacing();

        c_vector<double, 6> box_domain{};
        for (unsigned dim = 0; dim < 3; ++dim)
        {
            box_domain[2 * dim] = -2.0 * cell_radius;
            box_domain[2 * dim + 1] = 2.0 * cell_radius;
        }
        p_mesh->SetUpBoxCollection(interaction_cutoff, box_domain);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<3> cell_population(*p_mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), num_nodes);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1u);

        const double initial_gyration_radius = CalculateRadiusOfGyration(cell_population);
        TS_ASSERT(initial_gyration_radius > 0.0);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSemBasedSimulationWithSphericalElement");
        simulator.SetDt(0.005);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(0.05);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<3>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

        MAKE_PTR(SemForce<3>, p_sem_force);
        const double rho = 5.0;
        p_sem_force->SetIntraScalingFactor(rho);
        p_sem_force->SetInterScalingFactor(rho);
        const SemNScaledParameters intra_params = p_sem_force->ApplyNScaledIntraParameters(
            num_nodes, cell_radius, 20.0, 0.0, packing_density, 1.0);
        p_sem_force->SetIntraCutOffDistance(interaction_cutoff);
        p_sem_force->ApplyNScaledInterParameters(num_nodes, cell_radius, 20.0, 0.0, packing_density, 1.0);
        p_sem_force->SetInterCutOffDistance(interaction_cutoff);
        cell_population.SetDampingConstantNormal(intra_params.DampingConstant);
        simulator.AddForce(p_sem_force);

        /*
         * The scaled r_eq the force derives from N and the cell radius must match the spacing the
         * generator actually built, or the cell would start far from mechanical equilibrium.
         */
        TS_ASSERT_DELTA(p_sem_force->GetIntraEquilibriumDistance(), paper_r_eq, 1e-9);

        simulator.Solve();

        /*
         * The initial condition is a close packing, which is already at or near a minimum of the
         * pair potential, so the cell should barely move: surface nodes pull inward a little and
         * that is all. Assert it neither evaporates nor collapses. A wrong cut-off or a mismatched
         * r_eq shows up here as a large change in either direction.
         */
        const double final_gyration_radius = CalculateRadiusOfGyration(simulator.rGetCellPopulation());
        TS_ASSERT_DELTA(final_gyration_radius, initial_gyration_radius, 0.05 * initial_gyration_radius);

        // No node may leave the neighbourhood of the cell
        for (auto node_iter = simulator.rGetCellPopulation().rGetMesh().GetNodeIteratorBegin();
             node_iter != simulator.rGetCellPopulation().rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            TS_ASSERT_LESS_THAN(norm_2(node_iter->rGetLocation()), 2.0 * cell_radius);
        }
    }

    void TestSemBasedSimulationArchiving()
    {
        EXIT_IF_PARALLEL;  // Cell-based archiving does not work in parallel

        const unsigned num_nodes = 9u;  // 3x3 single element

        // Reference run: solve uninterrupted from 0 to 0.2 and record the final node positions.
        // A deterministic SemForce (no random force) is used so the archived run must match exactly.
        std::vector<c_vector<double, 2> > reference_positions(num_nodes);
        {
            SemSingleElementMeshGenerator<2> generator({3, 3}, 0.5);
            auto p_mesh = generator.GetMesh();
            c_vector<double, 4> domain{};
            domain[0] = -1.0; domain[1] = 2.0; domain[2] = -1.0; domain[3] = 2.0;
            p_mesh->SetUpBoxCollection(0.5, domain);

            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            SemBasedCellPopulation<2> cell_population(*p_mesh, cells);

            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestSemBasedSimulationArchiveReference");
            simulator.SetDt(0.01);
            simulator.SetSamplingTimestepMultiple(10);
            simulator.SetEndTime(0.2);
            simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2>>());
            simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

            // Use N-scaled parameters (packing=1 for the regular grid) so the Morse dynamics
            // are stable at the generator's node spacing.
            MAKE_PTR(SemForce<2>, p_sem_force);
            const SemNScaledParameters params = p_sem_force->ApplyNScaledIntraParameters(p_mesh->GetNumNodes(), 0.25, 20.0, 0.0, 1.0);
            p_sem_force->SetIntraCutOffDistance(0.5);
            cell_population.SetDampingConstantNormal(params.DampingConstant);
            simulator.AddForce(p_sem_force);

            simulator.Solve();

            for (unsigned i = 0; i < num_nodes; ++i)
            {
                reference_positions[i] = simulator.rGetCellPopulation().GetNode(i)->rGetLocation();
            }
        }

        // Reset the SimulationTime singleton so the save/restart run starts afresh at t = 0
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Save/restart run: solve from 0 to 0.1, checkpoint, then reload and continue to 0.2.
        // Saving exercises SemMesh's write-to-file archiving and reloading exercises the reader.
        {
            SemSingleElementMeshGenerator<2> generator({3, 3}, 0.5);
            auto p_mesh = generator.GetMesh();
            c_vector<double, 4> domain{};
            domain[0] = -1.0; domain[1] = 2.0; domain[2] = -1.0; domain[3] = 2.0;
            p_mesh->SetUpBoxCollection(0.5, domain);

            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            SemBasedCellPopulation<2> cell_population(*p_mesh, cells);

            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestSemBasedSimulationArchive");
            simulator.SetDt(0.01);
            simulator.SetSamplingTimestepMultiple(10);
            simulator.SetEndTime(0.1);
            simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2>>());
            simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

            MAKE_PTR(SemForce<2>, p_sem_force);
            const SemNScaledParameters params = p_sem_force->ApplyNScaledIntraParameters(p_mesh->GetNumNodes(), 0.25, 20.0, 0.0, 1.0);
            p_sem_force->SetIntraCutOffDistance(0.5);
            cell_population.SetDampingConstantNormal(params.DampingConstant);
            simulator.AddForce(p_sem_force);

            simulator.Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);
        }

        OffLatticeSimulation<2>* p_simulator =
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestSemBasedSimulationArchive", 0.1);

        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumNodes(), num_nodes);

        p_simulator->SetEndTime(0.2);
        p_simulator->Solve();

        // The restarted run must reproduce the uninterrupted reference run. The tolerance
        // allows for the six-significant-figure precision of the archived mesh node positions.
        for (unsigned i = 0; i < num_nodes; ++i)
        {
            c_vector<double, 2> loaded_position = p_simulator->rGetCellPopulation().GetNode(i)->rGetLocation();
            TS_ASSERT_DELTA(loaded_position[0], reference_positions[i][0], 1e-4);
            TS_ASSERT_DELTA(loaded_position[1], reference_positions[i][1], 1e-4);
        }

        delete p_simulator;
    }

};

#endif /*TESTSEMBASEDSIMULATION_HPP_*/
