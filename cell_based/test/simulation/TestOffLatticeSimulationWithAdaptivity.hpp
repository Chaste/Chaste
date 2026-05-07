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

#ifndef TESTOFFLATTICESIMULATIONWITHADAPTIVITY_HPP_
#define TESTOFFLATTICESIMULATIONWITHADAPTIVITY_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "RK4NumericalMethod.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Helper method to run a simple 3x3 NodeBased monolayer simulation starting from a compressed
 * hexagonal grid (natural spring length 1.0, initial grid spacing 0.1) so cells are
 * never at force equilibrium and the position trajectory is non-trivial.
 *
 * The population is set up with DifferentiatedCells (no divisions) so the
 * configuration is deterministic given the starting mesh and the numerical
 * method / time-step.
 *
 * Returns the final node positions sorted by initial index order.
 *
 * @param dt                     time step to use
 * @param useRk4                 if true use RK4NumericalMethod, else ForwardEuler
 * @param useAdaptive            if true enable adaptive time stepping
 * @param absoluteMovementThreshold  passed to SetAbsoluteMovementThreshold
 * @param outputDirectory        directory for VTK/results output
 */
static std::vector<c_vector<double,2> > RunMonolayer(
        double dt,
        bool useRk4,
        bool useAdaptive,
        double absoluteMovementThreshold,
        double endTime,
        const std::string& outputDirectory)
{
    // Ensure each simulation run starts from t=0 so output timestamps are comparable.
    SimulationTime::Destroy();
    SimulationTime::Instance()->SetStartTime(0.0);

    // -----------------------------------------------------------------------
    // Build a compressed 3x3 hexagonal grid (spacing 0.1, natural length 1.0)
    // so that cells experience net repulsive forces throughout the run.
    // -----------------------------------------------------------------------

    // Mesh
    const unsigned nx = 3;
    const unsigned ny = 3;
    const double mesh_spacing = 0.1;
    HoneycombMeshGenerator generator(nx, ny, 0);
    boost::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();
    p_generating_mesh->Scale(mesh_spacing, mesh_spacing);

    // Create NodesOnlyMesh - node-based populations don't remesh
    MAKE_PTR(NodesOnlyMesh<2>, p_mesh);
    p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

    // Differentiated cells – no divisions during the run
    std::vector<CellPtr> cells;
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
    CellsGenerator<FixedG1GenerationalCellCycleModel, 2> gen;
    gen.GenerateBasic(cells, p_mesh->GetNumNodes(), std::vector<unsigned>(), p_diff_type);

    // Population - NodeBasedCellPopulation never remeshes
    NodeBasedCellPopulation<2> population(*p_mesh, cells);
    population.SetAbsoluteMovementThreshold(absoluteMovementThreshold);

    // Simulation
    OffLatticeSimulation<2> sim(population);
    static unsigned run_counter = 0u;
    std::ostringstream output_dir;
    output_dir << outputDirectory << "_run_" << run_counter++;
    sim.SetOutputDirectory(output_dir.str());
    sim.SetDt(dt);
    sim.SetEndTime(endTime);
    sim.SetSamplingTimestepMultiple(1); // ensure results.viznodes is written at endTime

    // Force
    MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
    p_force->SetCutOffLength(1.5);
    p_force->SetMeinekeSpringStiffness(1.0);
    sim.AddForce(p_force);

    // Numerical method
    if (useRk4)
    {
        MAKE_PTR(RK4NumericalMethod<2>, p_method);
        p_method->SetUseAdaptiveTimestep(useAdaptive);
        sim.SetNumericalMethod(p_method);
    }
    else
    {
        MAKE_PTR(ForwardEulerNumericalMethod<2>, p_method);
        p_method->SetUseAdaptiveTimestep(useAdaptive);
        sim.SetNumericalMethod(p_method);
    }

    sim.Solve();

    // Collect final positions by node index for deterministic comparisons.
    std::vector<c_vector<double,2> > positions;
    positions.resize(population.GetNumNodes());
    for (auto node_iter = population.rGetMesh().GetNodeIteratorBegin();
         node_iter != population.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        positions[node_iter->GetIndex()] = node_iter->rGetLocation();
    }

    return positions;
}

/**
 * Compute the max-norm distance between two sets of node positions
 * (same ordering assumed).
 */
static double MaxPositionError(const std::vector<c_vector<double,2> >& a,
                               const std::vector<c_vector<double,2> >& b)
{
    assert(a.size() == b.size());
    double err = 0.0;
    for (unsigned i = 0; i < a.size(); ++i)
    {
        err = std::max(err, norm_2(a[i] - b[i]));
    }
    return err;
}

// ============================================================================

class TestOffLatticeSimulationWithAdaptivity : public AbstractCellBasedTestSuite
{
public:

    /**
     * TestFixedTimestepConvergence
     *
     * Runs the monolayer with ForwardEuler and RK4 at dt = 2^-4, ..., 2^-12
     * and compares to a reference RK4 solution at dt = 2^-14.  Checks that observed convergence orders are consistent with the expected first-order and fourth-order for FE and RK4 respectively, and that errors reduce as dt is reduced.
     */
    void TestFixedTimestepConvergence()
    {
        EXIT_IF_PARALLEL;

        const double threshold = 2.0; // default – no adaptive steps triggered
        const double common_dt = 1.0 / std::pow(2.0, 4.0);
        const double ref_dt    = 1.0 / std::pow(2.0, 14.0);

        // --------------------------------------------------------------------
        // Reference solution dt = 2^-14 and RK4 (should be very accurate)
        // --------------------------------------------------------------------
        const std::string ref_dir = "TestAdaptivity/Reference";
        std::vector<c_vector<double,2> > ref =
            RunMonolayer(ref_dt, /*rk4=*/true, /*adaptive=*/false, threshold, common_dt,
                         ref_dir);

        // --------------------------------------------------------------------
        // Collect errors at dt = 2^-4, ..., 2^-10
        // --------------------------------------------------------------------
        const unsigned n_levels = 6;
        std::vector<double> fe_errors(n_levels);
        std::vector<double> rk4_errors(n_levels);

        for (unsigned k = 0; k < n_levels; ++k)
        {
            double dt = 1.0 / std::pow(2.0, static_cast<double>(k + 4)); // 2^-4 ... 2^-10

            std::ostringstream fe_dir, rk4_dir;
            fe_dir  << "TestAdaptivity/FE_dt_2pow_neg"  << (k + 4);
            rk4_dir << "TestAdaptivity/RK4_dt_2pow_neg" << (k + 4);

            std::vector<c_vector<double,2> > fe_pos =
                RunMonolayer(dt, false, false, threshold, common_dt, fe_dir.str());
            std::vector<c_vector<double,2> > rk4_pos =
                RunMonolayer(dt, true,  false, threshold, common_dt, rk4_dir.str());

            fe_errors[k]  = MaxPositionError(fe_pos, ref);
            rk4_errors[k] = MaxPositionError(rk4_pos, ref);
        }

        // Data from a run with no adaptivity:
        // dt=2^-4 FE error=0.401622 RK4 error=0.0191423
        // dt=2^-5 FE error=0.112035 RK4 error=0.00134043
        // dt=2^-6 FE error=0.0427904 RK4 error=0.000130988
        // dt=2^-7 FE error=0.0190764 RK4 error=1.16722e-05
        // dt=2^-8 FE error=0.0090591 RK4 error=1.10516e-06
        // dt=2^-9 FE error=0.004421 RK4 error=8.18493e-08
        TS_ASSERT_DELTA(fe_errors[0], 0.401622, 1e-6);
        TS_ASSERT_DELTA(rk4_errors[0], 0.0191423, 1e-6);
        TS_ASSERT_DELTA(fe_errors[1], 0.112035, 1e-6);
        TS_ASSERT_DELTA(rk4_errors[1], 0.00134043, 1e-6);
        TS_ASSERT_DELTA(fe_errors[2], 0.0427904, 1e-6);
        TS_ASSERT_DELTA(rk4_errors[2], 0.000130988, 1e-6);
        TS_ASSERT_DELTA(fe_errors[3], 0.0190764, 1e-6);
        TS_ASSERT_DELTA(rk4_errors[3], 1.16722e-05, 1e-6);
        TS_ASSERT_DELTA(fe_errors[4], 0.0090591, 1e-6);
        TS_ASSERT_DELTA(rk4_errors[4], 1.10516e-06, 1e-6);
        TS_ASSERT_DELTA(fe_errors[5], 0.004421, 1e-6);
        TS_ASSERT_DELTA(rk4_errors[5], 8.18493e-08, 1e-6);


        // Report observed convergence order p from e(h)/e(h/2) = 2^p.
        // FE should be close to p=1 and RK4 close to p=4 on coarse refinements.
        double fe_order_8_to_9 = std::log(fe_errors[4] / fe_errors[5]) / std::log(2.0);
        double rk4_order_8_to_9 = std::log(rk4_errors[4] / rk4_errors[5]) / std::log(2.0);
        TS_ASSERT_DELTA(fe_order_8_to_9, 1.0, 0.1);
        TS_ASSERT_DELTA(rk4_order_8_to_9, 4.0, 0.5);


        // Some rough sanity checks on the errors:

        // RK4 errors should be smaller than FE at the same dt
        for (unsigned k = 0; k < n_levels; ++k)
        {
            TS_ASSERT_LESS_THAN(rk4_errors[k], fe_errors[k]);
        }

        // Each halving of dt should reduce the FE error by ~2× (first order) note we use 1.5 to account for numerical noise.
        for (unsigned k = 1; k < n_levels; ++k)
        {
            TS_ASSERT_LESS_THAN(1.5*fe_errors[k], fe_errors[k - 1]);
        }

        // Each halving of dt should reduce the RK4 error by ~16× (fourth order).
        // In practice numerical noise limits this so we are using 10.
        for (unsigned k = 1; k < n_levels; ++k)
        {

            TS_ASSERT_LESS_THAN(10*rk4_errors[k], rk4_errors[k - 1]);
        }

    }

    /**
     * TestAdaptiveForwardEuler
     *
     * Runs the monolayer with ForwardEuler and adaptive time-stepping starting
     * from a dt that triggers the Movement threshold.
     *
     */
    void TestAdaptiveForwardEuler()
    {
        EXIT_IF_PARALLEL;

        const double initial_dt = 1.0 / std::pow(2.0, 4.0);
        const double reference_amt = 2.0;  // Permissive for non-adaptive reference
        const std::vector<double> amts = {0.1, 0.01, 0.001};

        // Reference at fine dt (non-adaptive RK4) - use large threshold so it completes
        std::vector<c_vector<double,2> > ref =
            RunMonolayer(1.0 / std::pow(2.0, 14.0), /*rk4=*/true, /*adaptive=*/false,
                         /*absoluteMovementThreshold=*/reference_amt, initial_dt,
                         "TestAdaptivity/AdaptiveFERef");
        std::vector<double> errors(amts.size(), 0.0);
        for (unsigned i = 0; i < amts.size(); ++i)
        {
            std::ostringstream out_dir;
            out_dir << "TestAdaptivity/AdaptiveFE_amt_" << i;
            std::vector<c_vector<double,2> > adaptive =
                RunMonolayer(initial_dt, /*rk4=*/false, /*adaptive=*/true, amts[i], initial_dt,
                             out_dir.str());
            errors[i] = MaxPositionError(adaptive, ref);
        }

        // Data from a run with adaptive FE:
        // AMT=0.1  error=0.0190764
        // AMT=0.01 error=0.00108604
        // AMT=0.001 error=0.000135077
        TS_ASSERT_DELTA(errors[0], 0.0190764,   1e-6);
        TS_ASSERT_DELTA(errors[1], 0.00108604,  1e-6);
        TS_ASSERT_DELTA(errors[2], 0.000135077, 1e-6);
    }

    /**
     * TestAdaptiveRK4
     *
     * Same as TestAdaptiveForwardEuler but using the RK4 numerical method.
     * Because RK4 is more accurate per step, the error vs the reference
     * should be even smaller than for the adaptive FE case.
     */
    void TestAdaptiveRK4()
    {
        EXIT_IF_PARALLEL;

        const double initial_dt = 1.0 / std::pow(2.0, 4.0);
        const double reference_amt = 2.0;  // Permissive for non-adaptive reference
        const std::vector<double> amts = {0.1, 0.01, 0.001};

        // Reference at fine dt (non-adaptive RK4) - use large threshold so it completes
        std::vector<c_vector<double,2> > ref =
            RunMonolayer(1.0 / std::pow(2.0, 14.0), /*rk4=*/true, /*adaptive=*/false,
                         /*absoluteMovementThreshold=*/reference_amt, initial_dt,
                         "TestAdaptivity/AdaptiveRK4Ref");

        std::vector<double> errors(amts.size(), 0.0);
        for (unsigned i = 0; i < amts.size(); ++i)
        {
            std::ostringstream out_dir;
            out_dir << "TestAdaptivity/AdaptiveRK4_amt_" << i;

            std::vector<c_vector<double,2> > adaptive =
                RunMonolayer(initial_dt, /*rk4=*/true, /*adaptive=*/true, amts[i], initial_dt,
                             out_dir.str());

            errors[i] = MaxPositionError(adaptive, ref);
        }

        // Data from a run with adaptive RK4:
        // AMT=0.1  error=1.16722e-05
        // AMT=0.01 error=3.49773e-10
        // AMT=0.001 error=1.24127e-16
        TS_ASSERT_DELTA(errors[0], 1.16722e-05,  1e-9);
        TS_ASSERT_DELTA(errors[1], 3.49773e-10,  1e-12);
        TS_ASSERT_DELTA(errors[2], 1.24127e-16,  1e-16);// note this is close to machine precision for double, so we are using 1e-16 to allow for some numerical noise
    }
};

#endif // TESTOFFLATTICESIMULATIONWITHADAPTIVITY_HPP_
