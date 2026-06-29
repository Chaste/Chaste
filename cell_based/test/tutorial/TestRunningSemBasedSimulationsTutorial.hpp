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
/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTRUNNINGSEMBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGSEMBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * ## Examples showing how to create, run and visualize subcellular element method simulations
 *
 * ### Introduction
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize simulations
 * using the Subcellular Element Method (SEM). In the SEM, each biological cell is represented
 * by a collection of subcellular nodes (a `SemElement`) that interact through a modified Morse
 * potential. Full details of the mechanical model can be found in Sandersius & Newman, "Modelling
 * cell rheology with the Subcellular Element Model", Physical Biology, Vol. 5, No. 1, 2008.
 * doi:[10.1088/1478-3975/5/1/015002](https://doi.org/10.1088/1478-3975/5/1/015002).
 *
 * ### The test
 *
 * As in other cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

/* The following header is usually included in all cell-based test suites. */
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

/* We include the cell generators and cycle model required for all SEM simulations. SEM cells
 * do not divide by default, so we use `NoCellCycleModel`. */
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"

/* The simulation driver and smart pointer utilities. */
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "ForwardEulerNumericalMethod.hpp"

/* SEM mesh generators. `SemSingleElementMeshGenerator` creates a single cell from a regular
 * grid of nodes. `SemMultiElementMeshGenerator` creates a regular lattice of cells. */
#include "SemSingleElementMeshGenerator.hpp"
#include "SemMultiElementMeshGenerator.hpp"

/* The SEM cell population. */
#include "SemBasedCellPopulation.hpp"

/* SEM force laws. `SemForce` implements the modified Morse potential from Sandersius & Newman
 * (2008) with separate intra- and inter-cellular parameters. The two random force classes
 * provide the overdamped Langevin noise term: `SemGaussianRandomForce` generates independent
 * (uncorrelated) noise at each node, while `SemSpatiallyCorrelatedRandomForce` generates
 * noise that is spatially correlated between nearby nodes. */
#include "SemForce.hpp"
#include "SemGaussianRandomForce.hpp"
#include "SemSpatiallyCorrelatedRandomForce.hpp"

/* `SemParameterScaling.hpp` provides `SemComputeNScaledParameters<DIM>()`, which computes
 * a physically consistent parameter set (equilibrium distance, well depth, spring constant,
 * damping constant) from the number of nodes N, cell radius, and reference spring constant
 * κ₀, following the N-scaling laws in Sandersius & Newman (2008) Section 2. The two
 * convenience methods `ApplyNScaledIntraParameters()` and `ApplyNScaledInterParameters()`
 * on `SemForce` call this function and set the corresponding member variables directly. */
#include "SemParameterScaling.hpp"

/* Output writers that label each VTK point by the SEM element it belongs to and by its
 * region within that element. These are useful for visualizing cell boundaries in Paraview. */
#include "ElementIdNodePointDataWriter.hpp"
#include "NodeRegionPointDataWriter.hpp"

/*
 * Next, we define the test class.
 */
class TestRunningSemBasedSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    /*
     * ### Test 1 - a 2D single-cell SEM simulation with uncorrelated noise
     *
     * In the first test we run a simple SEM simulation of a single cell in 2D.
     * Each biological cell is represented by a 3×3 grid of subcellular nodes held
     * together by a modified Morse potential (`SemForce`). Thermal-like fluctuations
     * are added via `SemGaussianRandomForce`, which draws an independent Gaussian
     * random displacement for each node at every time step.
     */
    void TestSingleCell2D()
    {
        /* We begin by creating the SEM mesh. `SemSingleElementMeshGenerator<2>` constructs
         * a single `SemElement` containing a 3×3 regular grid of nodes (9 nodes total).
         * The second argument is the scale factor: the target node spacing in the
         * x-direction is `scaleFactor / numNodes[0]`, so here nodes are spaced
         * approximately 0.167 apart. */
        SemSingleElementMeshGenerator<2> generator({3, 3}, 0.5);
        auto p_mesh = generator.GetMesh();

        /* The SEM population uses a box collection to find interacting node pairs
         * efficiently. We must set up the box collection before constructing the
         * population. The first argument is the interaction cut-off distance (nodes
         * further apart than this are ignored by the force), and the second is the
         * bounding box [x_min, x_max, y_min, y_max] of the simulation domain. The
         * domain must be large enough to contain all nodes throughout the simulation. */
        const double interaction_cutoff = 0.25;
        c_vector<double, 4> box_domain;
        box_domain[0] = -1.0;
        box_domain[1] =  2.0;
        box_domain[2] = -1.0;
        box_domain[3] =  2.0;
        p_mesh->SetUpBoxCollection(interaction_cutoff, box_domain);

        /* Having created a mesh, we create one `CellPtr` per `SemElement`. Because SEM cells
         * do not divide we use `NoCellCycleModel`. */
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        /* We create the cell population, passing the mesh and cells. We also set the
         * node damping constant η used in the Langevin equation F = η * (velocity). */
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.0);

        /* We attach per-node data writers so that Paraview output includes the element
         * ID and region of each subcellular node. These are optional but make it easy
         * to colour-code cells during visualisation. */
        cell_population.AddNodePointDataWriter<ElementIdNodePointDataWriter>();
        cell_population.AddNodePointDataWriter<NodeRegionPointDataWriter>();

        /* We create an `OffLatticeSimulation` and set output directory, time step,
         * sampling interval, and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("SemSingleCell2D");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(1.0);

        /* SEM simulations require the Forward Euler numerical method. We must also
         * call `SetUseUpdateNodeLocation(false)` so that forces are applied via the
         * standard ODE integration pathway, not the alternative node-update path. */
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

        /* We add `SemForce` and configure its parameters using the N-dependent scaling
         * from Sandersius & Newman (2008) Section 2. `ApplyNScaledIntraParameters`
         * computes the equilibrium distance r_eq = 2 R_cell (p/N)^{1/DIM} and the Morse
         * well depth u₀ = κ r_eq² / (8 ρ²) from three physical inputs: the number of
         * subcellular nodes N, the cell radius R_cell, and a reference spring constant κ₀.
         * This ensures that if N is later changed to improve resolution, the macroscopic
         * cell stiffness is preserved automatically. For a single cell only the intra-cellular
         * parameters are relevant. The cut-off distance is set separately.
         *
         * The cell radius R_cell is half the scale factor (the half-width of the cell).
         * Passing packingDensity=1.0 tells the formula to treat the nodes as a regular
         * square grid — the correct packing for the generators used here. This choice
         * gives r_eq = scaleFactor / numNodes[0], which exactly matches the inter-node
         * spacing produced by SemSingleElementMeshGenerator. The reference spring constant
         * κ₀ controls the overall stiffness. */
        const double cell_radius = 0.25;   // half the scale factor
        const double kappa0 = 20.0;
        const double packing = 1.0;        // regular square grid

        MAKE_PTR(SemForce<2>, p_sem_force);
        p_sem_force->ApplyNScaledIntraParameters(p_mesh->GetNumNodes(), cell_radius, kappa0, 0.0, packing);
        p_sem_force->SetIntraCutOffDistance(interaction_cutoff);
        simulator.AddForce(p_sem_force);

        /* We add `SemGaussianRandomForce`, which applies independent Gaussian noise
         * to each node. The diffusion constant D controls the magnitude: the force
         * amplitude scales as η * sqrt(2D/dt). D must be small enough that the
         * resulting per-step displacements are much smaller than the node spacing,
         * otherwise the stiff Morse repulsion can overshoot and destabilise the
         * simulation. */
        MAKE_PTR(SemGaussianRandomForce<2>, p_noise);
        p_noise->SetDiffusionConstant(1e-6);
        simulator.AddForce(p_noise);

        /* To run the simulation, we call `Solve()`. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * If different simulation input parameters are being explored the lines should be
         * removed. */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1.0, 1e-10);
    }

    /*
     * To visualize the results, open Paraview and load
     * `$CHASTE_TEST_OUTPUT/SemSingleCell2D/results_from_time_0/results.pvd`.
     * Colour points by the `element_id` or `node_region` arrays to see the subcellular
     * node structure within the cell. See the
     * [Visualizing With Paraview](../visualizingwithparaview/) tutorial for more information.
     *
     * ### Test 2 - a 2D multi-cell SEM simulation with spatially correlated noise
     *
     * In the second test we simulate two SEM cells side by side in 2D. We now use
     * `SemMultiElementMeshGenerator` to place cells on a regular lattice and enable
     * inter-cellular interactions in `SemForce`. The noise model is changed to
     * `SemSpatiallyCorrelatedRandomForce`, which generates noise fields that are
     * correlated over a specified length scale. Spatially correlated noise causes nearby
     * nodes to fluctuate coherently, producing more realistic collective cell motion.
     */
    void TestMultiCell2D()
    {
        /* `SemMultiElementMeshGenerator<2>` creates a 2×1 lattice of cells (2 cells in
         * the x-direction, 1 in y). Each element contains a 3×3 grid of nodes, giving
         * 2 elements and 18 nodes in total. */
        SemMultiElementMeshGenerator<2> generator({3, 3}, {2, 1}, 0.5);
        auto p_mesh = generator.GetMesh();

        const double interaction_cutoff = 0.25;
        c_vector<double, 4> box_domain;
        box_domain[0] = -1.0;
        box_domain[1] =  2.0;
        box_domain[2] = -1.0;
        box_domain[3] =  2.0;
        p_mesh->SetUpBoxCollection(interaction_cutoff, box_domain);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.0);
        cell_population.AddNodePointDataWriter<ElementIdNodePointDataWriter>();
        cell_population.AddNodePointDataWriter<NodeRegionPointDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("SemMultiCell2D");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(1.0);

        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

        /* For a multi-cell simulation we set both intra-cellular (nodes within the same
         * element) and inter-cellular (nodes in different elements) force parameters.
         * N-scaling is applied per-cell: N is the number of nodes per cell, not the total
         * node count across the mesh. `ApplyNScaledInterParameters` uses the same scaling
         * formula as for intra-cellular interactions; in practice the inter-cellular κ₀
         * can be reduced independently to model weaker cell-cell adhesion. */
        const unsigned nodes_per_cell = p_mesh->GetNumNodes() / p_mesh->GetNumElements();
        const double cell_radius = 0.25;   // half the scale factor
        const double kappa0 = 20.0;
        const double packing = 1.0;        // regular square grid

        MAKE_PTR(SemForce<2>, p_sem_force);
        p_sem_force->ApplyNScaledIntraParameters(nodes_per_cell, cell_radius, kappa0, 0.0, packing);
        p_sem_force->SetIntraCutOffDistance(interaction_cutoff);
        p_sem_force->ApplyNScaledInterParameters(nodes_per_cell, cell_radius, kappa0, 0.0, packing);
        p_sem_force->SetInterCutOffDistance(interaction_cutoff);
        simulator.AddForce(p_sem_force);

        /* `SemSpatiallyCorrelatedRandomForce` requires the domain bounds and a correlation
         * length. Nodes closer together than `correlationLength` will receive similar noise
         * vectors at each step, causing coherent fluctuations within and between cells.
         * Setting `correlationLength` to the N-scaled equilibrium distance correlates
         * nearest neighbours; using a larger value produces longer-range coherence. */
        MAKE_PTR(SemSpatiallyCorrelatedRandomForce<2>, p_noise);
        p_noise->SetDiffusionConstant(1e-5);
        p_noise->SetCorrelationLength(p_sem_force->GetIntraEquilibriumDistance());
        p_noise->SetLowerCorner({{-1.0, -1.0}});
        p_noise->SetUpperCorner({{ 2.0,  2.0}});
        simulator.AddForce(p_noise);

        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial. */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 2u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1.0, 1e-10);
    }

    /*
     * To visualize the results, load
     * `$CHASTE_TEST_OUTPUT/SemMultiCell2D/results_from_time_0/results.pvd`
     * in Paraview and colour points by `element_id` to distinguish the two cells.
     *
     * ### Test 3 - a 3D single-cell SEM simulation with spatially correlated noise
     *
     * In the third test we run a single-cell SEM simulation in 3D. The setup is very
     * similar to the 2D test, with the dimension template changed from 2 to 3 and a
     * 3D bounding box. We continue to use `SemSpatiallyCorrelatedRandomForce` to
     * illustrate how the same correlated noise model extends naturally to 3D.
     */
    void TestSingleCell3D()
    {
        /* `SemSingleElementMeshGenerator<3>` creates a 3×3×3 grid of nodes (27 nodes)
         * forming a single 3D SEM element. */
        SemSingleElementMeshGenerator<3> generator({3, 3, 3}, 0.5);
        auto p_mesh = generator.GetMesh();

        /* In 3D the box collection domain vector has six entries:
         * [x_min, x_max, y_min, y_max, z_min, z_max]. */
        const double interaction_cutoff = 0.25;
        c_vector<double, 6> box_domain;
        box_domain[0] = -1.0;
        box_domain[1] =  2.0;
        box_domain[2] = -1.0;
        box_domain[3] =  2.0;
        box_domain[4] = -1.0;
        box_domain[5] =  2.0;
        p_mesh->SetUpBoxCollection(interaction_cutoff, box_domain);

        /* We compute the full N-scaled parameter set using `SemComputeNScaledParameters`.
         * This returns all four N-dependent quantities at once: the equilibrium distance
         * r_eq, the Morse well depth u₀, the harmonic spring constant κ, and the damping
         * constant η = η₀ N. Setting the cell population's damping constant to the N-scaled
         * value η ensures that the Langevin equation η ẋ = F_det + ξ is correctly balanced
         * as N changes. The rho parameter must match the force's scaling factor (default 5.0);
         * η₀ is the per-node reference damping. Here we choose η₀ = 1/N so that the
         * total damping η = η₀ N = 1 is independent of N; different N values can then be
         * compared directly without re-tuning the diffusion constant. */
        const unsigned num_nodes = p_mesh->GetNumNodes();
        const double cell_radius = 0.25;   // half the scale factor
        const double kappa0 = 20.0;
        const double rho = 5.0;
        const double packing = 1.0;        // regular cubic grid
        const double eta0 = 1.0 / static_cast<double>(num_nodes);
        const SemNScaledParameters nscaled =
            SemComputeNScaledParameters<3>(num_nodes, cell_radius, kappa0, rho, 0.0, eta0, packing);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        SemBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(nscaled.DampingConstant);
        cell_population.AddNodePointDataWriter<ElementIdNodePointDataWriter>();
        cell_population.AddNodePointDataWriter<NodeRegionPointDataWriter>();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("SemSingleCell3D");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(1.0);

        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<3>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

        /* We configure the force using the same N-scaling parameters. Setting
         * `IntraScalingFactor` to `rho` before calling `ApplyNScaledIntraParameters`
         * ensures both the force and the damping calculation use the same ρ. For a
         * single cell only intra-cellular interactions are relevant. */
        MAKE_PTR(SemForce<3>, p_sem_force);
        p_sem_force->SetIntraScalingFactor(rho);
        p_sem_force->ApplyNScaledIntraParameters(num_nodes, cell_radius, kappa0, 0.0, packing);
        p_sem_force->SetIntraCutOffDistance(interaction_cutoff);
        simulator.AddForce(p_sem_force);

        /* For the 3D spatially correlated noise we must supply 3D corners. The
         * correlation length is set to the N-scaled equilibrium distance, so coherent
         * fluctuations couple nearest-neighbour nodes. */
        MAKE_PTR(SemSpatiallyCorrelatedRandomForce<3>, p_noise);
        p_noise->SetDiffusionConstant(1e-5);
        p_noise->SetCorrelationLength(nscaled.EquilibriumDistance);
        p_noise->SetLowerCorner({{-1.0, -1.0, -1.0}});
        p_noise->SetUpperCorner({{ 2.0,  2.0,  2.0}});
        simulator.AddForce(p_noise);

        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial. */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1.0, 1e-10);
    }
};
/*
 * To visualize the 3D results, use Paraview. See the
 * [Visualizing With Paraview](../visualizingwithparaview/) tutorial for more information.
 *
 * Load the file `$CHASTE_TEST_OUTPUT/SemSingleCell3D/results_from_time_0/results.pvd`,
 * add a Glyph filter with Sphere representation, and colour by `element_id` or
 * `node_region` to see the subcellular structure of the cell.
 */

#endif /* TESTRUNNINGSEMBASEDSIMULATIONSTUTORIAL_HPP_ */
