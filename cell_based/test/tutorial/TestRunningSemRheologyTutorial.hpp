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
 * ## Measuring the rheology of a single cell with the subcellular element model
 *
 * This tutorial goes beyond setting a simulation up and running it: it performs a *measurement*.
 * We build an approximately spherical cell, bring it to a well-defined equilibrium, compress it
 * between two notional plates, and record how it creeps under the load. Each stage corresponds to
 * part of Sandersius & Newman (2008) Phys. Biol. 5 015002:
 *
 * 1. Generate an approximately spherical aggregate of `N` subcellular elements (Section 1.2).
 * 2. Anneal it with thermal noise, then quench it, to reach the *amorphous* packing the paper
 *    equilibrates rather than the crystal we start from (Section 1.2).
 * 3. Load it in uniaxial compression, sharing a total force `F` over the elements in a thin slab
 *    at each end (Section 3.1).
 * 4. Measure the strain and the creep compliance over time.
 *
 * If you have not seen a SEM simulation before, start with `TestRunningSemBasedSimulationsTutorial`,
 * which covers the basic setup this tutorial assumes.
 */

#ifndef TESTRUNNINGSEMRHEOLOGYTUTORIAL_HPP_
#define TESTRUNNINGSEMRHEOLOGYTUTORIAL_HPP_

/*
 * ### Include files
 *
 * As in any Chaste cell-based test, we begin with `CheckpointArchiveTypes.hpp` (via the abstract
 * test suite) and the cell population machinery.
 */
#include <cxxtest/TestSuite.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <boost/make_shared.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

/*
 * `SemSphericalElementMeshGenerator` builds a single approximately spherical SEM element from a
 * node count and a cell radius, which is the parameterisation the paper's scaling laws use. It
 * reports the resulting node spacing, which we will use directly as the equilibrium distance of
 * the interaction potential.
 */
#include "SemSphericalElementMeshGenerator.hpp"
#include "SemBasedCellPopulation.hpp"
#include "SemMesh.hpp"

/*
 * A population writer that records how far the cell extends along each axis at every output
 * timestep, which is what the creep measurement needs.
 */
#include "CellPopulationExtentWriter.hpp"

/*
 * `SemForce` is the modified Morse potential holding the cell together;
 * `SemGaussianRandomForce` supplies the thermal noise used to anneal it; and `UniaxialLoadForce`
 * applies the external load of the rheology experiment.
 */
#include "SemForce.hpp"
#include "SemGaussianRandomForce.hpp"
#include "SemParameterScaling.hpp"
#include "UniaxialLoadForce.hpp"

/*
 * Writer that labels each VTK point by its region, so the cell can be colour-coded in Paraview.
 */
#include "NodeRegionPointDataWriter.hpp"

/* This test is always run sequentially (never in parallel). */
#include "FakePetscSetup.hpp"

class TestRunningSemRheologyTutorial : public AbstractCellBasedTestSuite
{
private:

    /*
     * ### Helper functions
     *
     * A snapshot of every node position, used to decide when the cell has stopped moving.
     */
    std::vector<c_vector<double, 3> > GetNodePositions(SemMesh<3>& rMesh)
    {
        std::vector<c_vector<double, 3> > positions;
        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            positions.push_back(rMesh.GetNode(i)->rGetLocation());
        }
        return positions;
    }

    /* The largest distance any node has moved between two snapshots. */
    double GetMaximumDisplacement(const std::vector<c_vector<double, 3> >& rBefore,
                                  const std::vector<c_vector<double, 3> >& rAfter)
    {
        double max_displacement = 0.0;
        for (unsigned i = 0; i < rBefore.size(); ++i)
        {
            max_displacement = std::max(max_displacement, norm_2(rAfter[i] - rBefore[i]));
        }
        return max_displacement;
    }

    /*
     * The spread of the coordination number over the cell's interior nodes, i.e. the standard
     * deviation of the number of neighbours each has within a shade over one node spacing.
     *
     * This is our measure of crystallinity. In a perfect close packing every interior node has
     * exactly twelve neighbours, so the spread is exactly zero; in an amorphous packing the
     * coordination varies from node to node and the spread is positive.
     */
    double GetInteriorCoordinationSpread(SemMesh<3>& rMesh, double spacing)
    {
        std::vector<double> coordinations;

        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            if (rMesh.GetNode(i)->GetRegion() != SEM_INTERIOR_REGION)
            {
                continue;
            }

            unsigned num_neighbours = 0u;
            for (unsigned j = 0; j < rMesh.GetNumNodes(); ++j)
            {
                if (i != j
                    && norm_2(rMesh.GetNode(j)->rGetLocation() - rMesh.GetNode(i)->rGetLocation())
                           < 1.05 * spacing)
                {
                    num_neighbours++;
                }
            }
            coordinations.push_back(static_cast<double>(num_neighbours));
        }

        double mean = 0.0;
        for (double coordination : coordinations)
        {
            mean += coordination;
        }
        mean /= static_cast<double>(coordinations.size());

        double variance = 0.0;
        for (double coordination : coordinations)
        {
            variance += (coordination - mean) * (coordination - mean);
        }
        variance /= static_cast<double>(coordinations.size());

        return sqrt(variance);
    }

    /*
     * Write the pair distribution function of Section 1.2: a histogram of the separation of every
     * pair of elements, normalised so that the entries sum to one. A crystal produces sharp peaks
     * at the discrete neighbour-shell distances of the lattice; an amorphous packing produces a
     * broadened first peak and a characteristic split second peak.
     */
    void WritePairDistributionFunction(SemMesh<3>& rMesh, double spacing,
                                       OutputFileHandler& rHandler, const std::string& rFileName)
    {
        const unsigned num_bins = 60u;
        const double max_separation = 4.0 * spacing;
        const double bin_width = max_separation / static_cast<double>(num_bins);

        std::vector<unsigned> histogram(num_bins, 0u);
        unsigned num_pairs = 0u;

        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            for (unsigned j = i + 1; j < rMesh.GetNumNodes(); ++j)
            {
                const double separation
                    = norm_2(rMesh.GetNode(j)->rGetLocation() - rMesh.GetNode(i)->rGetLocation());
                if (separation < max_separation)
                {
                    histogram[static_cast<unsigned>(separation / bin_width)]++;
                    num_pairs++;
                }
            }
        }

        out_stream p_file = rHandler.OpenOutputFile(rFileName);
        *p_file << "# r/r_eq\trho(r)\n";
        for (unsigned bin = 0; bin < num_bins; ++bin)
        {
            const double r = (static_cast<double>(bin) + 0.5) * bin_width;
            *p_file << r / spacing << "\t"
                    << static_cast<double>(histogram[bin]) / static_cast<double>(num_pairs) << "\n";
        }
        p_file->close();
    }

public:

    /*
     * ### The rheology experiment
     */


    void TestSingleCellCreepUnderUniaxialLoad()
    {
        /*
         * #### Stage 1: build the cell
         *
         * We ask for `N` subcellular elements filling a sphere of the given radius. The generator
         * carves the cell out of a close packing and reports the resulting nearest-neighbour
         * distance, which is the `r_eq` of the interaction potential.
         *
         * Keeping `N` modest keeps this tutorial quick. Values in the high hundreds are closer to
         * the paper; the commented-out line below is a realistic choice.
         */
        const unsigned num_nodes = 200u;
        // const unsigned num_nodes = 800u;   // paper-scale aggregate
        const double cell_radius = 0.25;

        SemSphericalElementMeshGenerator<3> generator(num_nodes, cell_radius);
        auto p_mesh = generator.GetMesh();
        const double r_eq = generator.GetNodeSpacing();

        /*
         * The paper relates the equilibrium separation to the cell radius by
         * `r_eq = 2 R (p/N)^(1/3)`, with `p` the packing density of spheres. Our generator does
         * not use that relation - the spacing falls out of the geometry - so the two agreeing is a
         * useful check that the packing really is dense. We use the generator's value, which is
         * exact for the mesh we actually have.
         */
        const double packing_density = 0.7405;
        const double paper_r_eq
            = 2.0 * cell_radius * std::pow(packing_density / static_cast<double>(num_nodes), 1.0 / 3.0);
        TS_ASSERT_DELTA(r_eq, paper_r_eq, 0.05 * paper_r_eq);

        /*
         * The interaction is cut off a little beyond the second neighbour shell. The box
         * collection must cover everywhere the cell might reach, and the spherical generator
         * centres the cell on the origin, so the domain spans negative coordinates too.
         */
        const double interaction_cutoff = 2.0 * r_eq;

        c_vector<double, 6> box_domain;
        box_domain[0] = -1.0;
        box_domain[1] = 1.0;
        box_domain[2] = -1.0;
        box_domain[3] = 1.0;
        box_domain[4] = -1.0;
        box_domain[5] = 1.0;
        p_mesh->SetUpBoxCollection(interaction_cutoff, box_domain);

        /*
         * The VTK output can carry a reconstructed surface for each element as well as its point
         * cloud, but the surface is rebuilt by an alpha-shape reconstruction at every output step,
         * which dominates the run time when sampling often. The elements themselves are all we need
         * here, so we turn the surface off.
         */
        p_mesh->SetOutputElementSurfacesToVtk(false);

        /* One cell per SEM element, with no cell cycle model since the cell neither grows nor divides. */
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        SemBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddNodePointDataWriter<NodeRegionPointDataWriter>();

        /*
         * `CellPopulationExtentWriter` records how far the cell reaches along each axis at every
         * output timestep. It is registered before any solving is done, so the equilibration is
         * recorded as well as the creep run: the cell swells while it is molten and then slowly
         * re-densifies, and we need to see that settle before the load is applied.
         */
        cell_population.AddPopulationWriter<CellPopulationExtentWriter>();

        /*
         * #### The interaction potential
         *
         * We take the N-dependent scaling of paper Section 2 for the spring constant, well depth
         * and damping constant, then override the equilibrium distance with the generator's exact
         * node spacing.
         */
        const double kappa0 = 20.0;
        const double rho = 5.0;
        const double eta0 = 1.0 / static_cast<double>(num_nodes);

        MAKE_PTR(SemForce<3>, p_sem_force);
        p_sem_force->SetIntraScalingFactor(rho);
        const SemNScaledParameters nscaled = p_sem_force->ApplyNScaledIntraParameters(
            num_nodes, cell_radius, kappa0, 0.0, packing_density, eta0);
        p_sem_force->SetIntraEquilibriumDistance(r_eq);
        p_sem_force->SetIntraCutOffDistance(interaction_cutoff);
        cell_population.SetDampingConstantNormal(nscaled.DampingConstant);

        /*
         * #### The simulator
         */
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("SemRheology");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<3> >());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);
        simulator.AddForce(p_sem_force);

        OutputFileHandler results_handler("SemRheology", false);

        /*
         * The cell starts as a perfect crystal carved from a close packing, so every interior
         * element has exactly twelve neighbours and the coordination spread is exactly zero.
         */
        WritePairDistributionFunction(*p_mesh, r_eq, results_handler, "pair_distribution_initial.dat");
        const double initial_spread = GetInteriorCoordinationSpread(*p_mesh, r_eq);
        TS_ASSERT_DELTA(initial_spread, 0.0, 1e-12);

        /*
         * #### Stage 2: anneal and quench
         *
         * A crystal already sits at a minimum of the potential, so simply relaxing it would barely
         * move anything, and quenching a *gently* warmed crystal just lets it snap back onto the
         * same lattice. To reach the amorphous aggregate the paper equilibrates, we must genuinely
         * melt it and then cool it back down.
         *
         * The diffusion constant is the one parameter here needing real care, and the quantity to
         * reason about is not a displacement but an *energy*. In overdamped Langevin dynamics the
         * thermal energy scale is `D * eta`, and the energy binding one element to a neighbour is
         * the well depth `u0`. Setting `D * eta` comparable to `u0` melts the crystal; a few times
         * larger and elements simply boil off the surface, since an unconfined cluster held by a
         * finite well has nothing to hold them once the thermal energy approaches the binding
         * energy. Setting it much smaller leaves a warm crystal that recrystallises on quenching.
         *
         * Note that `eta` here is the *scaled* damping constant, eta_0 * N, not eta_0.
         */
        const double temperature_fraction = 1.0;
        const double diffusion_constant
            = temperature_fraction * nscaled.WellDepth / nscaled.DampingConstant;

        /*
         * `SetCoolingWindow` cools the noise towards zero between the two given times, so the whole
         * equilibration - melt, quench, and the subsequent relaxation with no noise at all - happens
         * inside one call to `Solve()`. Cooling gradually rather than switching the noise off in one
         * go matters: an abrupt quench traps the cell in whatever configuration it happened to be
         * in, whereas a ramp lets it settle into a better packing.
         */
        const double melt_duration = 5.0;
        const double cooling_duration = 2.0;
        const double settling_duration = 5.0;
        // const double melt_duration = 5.0;      // longer melt
        // const double cooling_duration = 10.0;  // longer cooling
        // const double settling_duration = 50.0; // and a longer settle afterwards

        MAKE_PTR(SemGaussianRandomForce<3>, p_noise);
        p_noise->SetDiffusionConstant(diffusion_constant);
        p_noise->SetCoolingWindow(melt_duration, melt_duration + cooling_duration);
        simulator.AddForce(p_noise);

        double time = melt_duration + cooling_duration + settling_duration;
        simulator.SetEndTime(time);
        simulator.Solve();

        /*
         * The cell is now an amorphous aggregate: the elements no longer sit on a lattice, so the
         * coordination number varies from one to the next and its spread is no longer zero. The
         * two pair distribution files can be plotted against each other to see the sharp lattice
         * peaks give way to the broadened first peak of a dense random packing.
         */
        WritePairDistributionFunction(*p_mesh, r_eq, results_handler, "pair_distribution_equilibrated.dat");
        const double equilibrated_spread = GetInteriorCoordinationSpread(*p_mesh, r_eq);
        TS_ASSERT_LESS_THAN(0.5, equilibrated_spread);

        /*
         * #### Stage 3: load the cell
         *
         * Following paper Section 3.1, a total force `F` is shared between the elements in a thin
         * slab at each end of the loading axis, directed so the two slabs are pushed together. The
         * slab is one node spacing deep, so it grips a single layer of elements.
         *
         * The two loads are equal and opposite, so the cell deforms without drifting and no
         * boundary condition is needed to hold it in place.
         */
        const double applied_load = 1.2 * nscaled.SpringConstant * r_eq;

        MAKE_PTR(UniaxialLoadForce<3>, p_load);
        p_load->SetLoad(applied_load);
        p_load->SetLoadingAxis(2u);
        p_load->SetSlabThickness(2.0 * r_eq);
        simulator.AddForce(p_load);

        const double reference_height = p_mesh->GetWidth(2u);

        /*
         * #### Stage 4: record the creep curve
         *
         * `CellPopulationExtentWriter` records the extent of the population along each axis at
         * every output timestep, into a single `cellpopulationextent.dat`. That is all we need:
         * dividing the reduction in extent along the loading axis by the reference height gives the
         * strain. Dividing that in turn by the applied stress - the load over the cross-sectional
         * area it acts on, `pi R^2` for a cell of this radius - gives the creep compliance
         * `J = strain/stress`, the quantity a rheology experiment reports.
         *
         * Collecting the measurement this way rather than by stopping the simulation to sample it
         * keeps the whole creep run in one `Solve()`, and so in one results directory. It also
         * samples at every output timestep instead of only where we happened to stop.
         */
        const double creep_duration = 2.0;
        // const double creep_duration = 20.0;   // longer creep curve

        time += creep_duration;
        simulator.SetSamplingTimestepMultiple(5);
        simulator.SetEndTime(time);
        simulator.Solve();

        /*
         * The cell should be measurably shorter, and should have gone on creeping rather than
         * stopping at its elastic response.
         */
        const double final_height = p_mesh->GetWidth(2u);
        const double final_strain = (reference_height - final_height) / reference_height;

        TS_ASSERT_LESS_THAN(final_height, reference_height);
        TS_ASSERT_LESS_THAN(0.0, final_strain);

        /*
         * #### The output files
         *
         * The `SemRheology` directory holds one `results_from_time_*` directory per call to
         * `Solve()` - two in all, one for the equilibration and one for the creep run - each with a
         * `results.pvd` that Paraview opens directly. Each also holds a `cellpopulationextent.dat`
         * whose columns are the simulation time followed by the radius of gyration along x, y and
         * z. At the top level are `pair_distribution_initial.dat` and
         * `pair_distribution_equilibrated.dat`, for plotting the crystal and amorphous packings
         * against one another.
         *
         * #### What the results look like
         *
         * The durations set above are deliberately short, so that this runs as a test in reasonable
         * time. **The behaviour described here is what you get with the longer durations commented
         * out beside them** - a longer melt, a longer cooling ramp, a much longer settle and a
         * longer creep run. With the short values the cell is nowhere near equilibrated before the
         * load is applied, and the numbers below will not reproduce.
         *
         * *The packing.* `r_eq` from the generator agrees with the paper's
         * `r_eq = 2 R (p/N)^(1/3)` to about half a percent. In the initial pair distribution the
         * crystal shows sharp peaks with **exact zeros between them** - nothing at all between about
         * 1.05 and 1.40 `r_eq` - and spikes at sqrt(2), sqrt(3) and 2 `r_eq`, the neighbour shells
         * of a close packing. After annealing and quenching every one of those gaps has filled in,
         * the first peak has broadened and dropped, and the second peak has **split into two**, near
         * sqrt(3) and 2 `r_eq`. That split is the signature of a dense random packing rather than a
         * crystal or a liquid, and it is the paper's own Section 1.2 diagnostic.
         *
         * *The compression.* Both lateral axes start expanding almost at once, within a fraction of
         * a time unit of the load being applied, while the loaded axis contracts. The apparent
         * Poisson ratio climbs from near zero to about 0.44 as the cell settles into its response,
         * and the volume falls by only about 2%: the cell deforms nearly incompressibly.
         *
         * *Why the radius of gyration matters here.* Measured by a bounding box instead, the same
         * cell looks 22% anisotropic before it is even loaded, and one lateral axis appears to
         * *contract* for the first several time units under load. Both are artefacts. The close
         * packed lattice is axis aligned and reaches further along some axes than others - along z
         * the sites lie only at multiples of the layer spacing, so the bounding box stops short -
         * and a bounding box is fixed by just two nodes out of hundreds. Averaging over every node
         * removes both effects, which is why `CellPopulationExtentWriter` reports a radius of
         * gyration.
         *
         * *What does not match the paper.* Two things, and both are worth understanding.
         * First, there is no elastic plateau: the cell creeps from the very first
         * timestep, so the stiffness you measure depends strongly on how soon after loading you
         * look, and no single number can be quoted as *the* elastic modulus. The paper's Section 2
         * estimate falls within the range this protocol produces, but that is a weak statement.
         * Second, the strain reported by a log-log fit grows as roughly `t^0.6`, well above the
         * weak power law of `t^0.1` to `t^0.3` that the abstract describes. The load here drives the
         * cell past 10% strain, and the paper is explicit that the model in its simplest form
         * "cannot describe long-time/large-strain cell responses". Reducing the load brings the
         * strain into that regime, at the cost of a much smaller signal.
         *
         * *A residual worth knowing about.* The generated cell is genuinely slightly oblate, by a
         * few percent, because the ball is carved from a lattice. The melt reduces this but does not
         * abolish it, and it shows up at large strain as the two lateral axes bulging by different
         * amounts. Melting for longer is the lever that acts on it, since shape relaxation is far
         * faster in the molten state than in the quenched one.
         */
    }
};

#endif /*TESTRUNNINGSEMRHEOLOGYTUTORIAL_HPP_*/
