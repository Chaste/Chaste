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

#ifndef TESTSEMPARAMETERSCALING_HPP_
#define TESTSEMPARAMETERSCALING_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>

#include "SemParameterScaling.hpp"
#include "SemForce.hpp"
#include "SemLinearForce.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Tests for SemComputeNScaledParameters and the ApplyNScaled* convenience methods
 * on SemForce and SemLinearForce.
 *
 * Reference: Sandersius & Newman (2008) Phys. Biol. 5 015002, Section 2.
 */
class TestSemParameterScaling : public CxxTest::TestSuite
{
public:

    /**
     * Verify SemComputeNScaledParameters against hand-calculated values for a
     * clean-integer case in 3D (N=27, so N^{1/3}=3 exactly).
     */
    void TestComputeNScaledParameters3DCleanInteger()
    {
        // N=27, DIM=3: N^{1/3}=3 exactly, n_factor = 1/3.
        // kappa0=9, lambda=1: kappa = 9 * (1/3) * (1 - 1/3) = 9*(1/3)*(2/3) = 2
        const unsigned N = 27u;
        const double rho = 5.0;
        const double r_cell = 1.0;
        const double kappa0 = 9.0;
        const double lambda = 1.0;
        const double eta0 = 2.0;
        const double packing = 0.7405;  // explicit FCC

        SemNScaledParameters p = SemComputeNScaledParameters<3>(N, r_cell, kappa0, rho, lambda, eta0, packing);

        const double r_eq_expected = 2.0 * std::pow(packing / 27.0, 1.0/3.0);
        TS_ASSERT_DELTA(p.EquilibriumDistance, r_eq_expected, 1e-9);

        TS_ASSERT_DELTA(p.SpringConstant, 2.0, 1e-9);

        const double u0_expected = 2.0 * r_eq_expected * r_eq_expected / (8.0 * 25.0);
        TS_ASSERT_DELTA(p.WellDepth, u0_expected, 1e-9);

        // Harmonic self-consistency: kappa = 8*rho^2*u0/r_eq^2
        TS_ASSERT_DELTA(8.0 * rho * rho * p.WellDepth / (p.EquilibriumDistance * p.EquilibriumDistance),
                        p.SpringConstant, 1e-9);

        TS_ASSERT_DELTA(p.DampingConstant, eta0 * N, 1e-9);
    }

    /**
     * Verify SemComputeNScaledParameters for lambda=0 in 3D with random close-packing
     * density (the value used in the Sandersius paper).
     */
    void TestComputeNScaledParameters3DLambdaZero()
    {
        const unsigned N = 100u;
        const double rho = 4.0;
        const double r_cell = 5.0;
        const double kappa0 = 1.0;
        const double packing = 0.64;  // random close packing (Sandersius 2008)

        SemNScaledParameters p = SemComputeNScaledParameters<3>(N, r_cell, kappa0, rho, 0.0, 1.0, packing);

        const double r_eq_expected = 2.0 * r_cell * std::pow(packing / N, 1.0/3.0);
        const double kappa_expected = kappa0 * std::pow(static_cast<double>(N), -1.0/3.0);
        const double u0_expected = kappa_expected * r_eq_expected * r_eq_expected / (8.0 * rho * rho);

        TS_ASSERT_DELTA(p.EquilibriumDistance, r_eq_expected, 1e-9);
        TS_ASSERT_DELTA(p.SpringConstant, kappa_expected, 1e-9);
        TS_ASSERT_DELTA(p.WellDepth, u0_expected, 1e-9);
        TS_ASSERT_DELTA(p.DampingConstant, static_cast<double>(N), 1e-9);

        TS_ASSERT_DELTA(8.0 * rho * rho * p.WellDepth / (p.EquilibriumDistance * p.EquilibriumDistance),
                        p.SpringConstant, 1e-9);
    }

    /**
     * Verify the 2D default packing density (HCP, π/(2√3) ≈ 0.9069).
     */
    void TestComputeNScaledParameters2DDefaultPacking()
    {
        const unsigned N = 16u;
        const double rho = 3.0;
        const double r_cell = 2.0;
        const double kappa0 = 1.0;
        const double packing_2d = 0.9069;  // HCP default

        SemNScaledParameters p = SemComputeNScaledParameters<2>(N, r_cell, kappa0, rho);

        const double r_eq_expected = 2.0 * r_cell * std::pow(packing_2d / N, 0.5);
        TS_ASSERT_DELTA(p.EquilibriumDistance, r_eq_expected, 1e-4);

        TS_ASSERT_DELTA(8.0 * rho * rho * p.WellDepth / (p.EquilibriumDistance * p.EquilibriumDistance),
                        p.SpringConstant, 1e-9);
    }

    /**
     * Verify that packing=1.0 with R_cell=scaleFactor/2 gives r_eq equal to the
     * inter-node spacing scaleFactor/numNodes produced by the regular-grid generators.
     * This is the canonical parameter choice for SemSingleElementMeshGenerator /
     * SemMultiElementMeshGenerator with a regular cubic/square lattice.
     */
    void TestRegularGridPackingGivesCorrectNodeSpacing()
    {
        const double scale_factor = 0.5;
        const double r_cell = scale_factor / 2.0;  // = 0.25
        const double packing = 1.0;                // regular grid

        // 2D: 3x3 = 9 nodes, spacing = scale_factor / 3 = 0.1667
        {
            const unsigned N = 9u;
            SemNScaledParameters p = SemComputeNScaledParameters<2>(N, r_cell, 1.0, 5.0, 0.0, 1.0, packing);
            const double expected_spacing = scale_factor / 3.0;
            TS_ASSERT_DELTA(p.EquilibriumDistance, expected_spacing, 1e-9);
        }

        // 3D: 3x3x3 = 27 nodes, spacing = scale_factor / 3 = 0.1667
        {
            const unsigned N = 27u;
            SemNScaledParameters p = SemComputeNScaledParameters<3>(N, r_cell, 1.0, 5.0, 0.0, 1.0, packing);
            const double expected_spacing = scale_factor / 3.0;
            TS_ASSERT_DELTA(p.EquilibriumDistance, expected_spacing, 1e-9);
        }
    }

    /**
     * Verify ApplyNScaledIntraParameters overwrites mIntraEquilibriumDistance and
     * mIntraWellDepth on SemForce, leaving mIntraScalingFactor (rho) unchanged.
     */
    void TestApplyNScaledIntraParametersOnSemForce()
    {
        SemForce<3> force;
        force.SetIntraScalingFactor(5.0);
        force.SetIntraWellDepth(99.0);           // sentinel — should be overwritten
        force.SetIntraEquilibriumDistance(99.0); // sentinel — should be overwritten

        // eta0 = 2.0, so the returned damping constant should be eta0*N = 54
        const SemNScaledParameters returned = force.ApplyNScaledIntraParameters(27u, 1.0, 9.0, 1.0, 0.7405, 2.0);

        // kappa = 9*(1/3)*(1 - 1/3) = 2
        const double r_eq = 2.0 * std::pow(0.7405 / 27.0, 1.0/3.0);
        const double u0 = 2.0 * r_eq * r_eq / (8.0 * 25.0);
        TS_ASSERT_DELTA(force.GetIntraEquilibriumDistance(), r_eq, 1e-9);
        TS_ASSERT_DELTA(force.GetIntraWellDepth(), u0, 1e-9);
        TS_ASSERT_DELTA(force.GetIntraScalingFactor(), 5.0, 1e-9);  // rho unchanged

        // The returned parameter set should match the force state and carry the damping constant
        TS_ASSERT_DELTA(returned.EquilibriumDistance, r_eq, 1e-9);
        TS_ASSERT_DELTA(returned.WellDepth, u0, 1e-9);
        TS_ASSERT_DELTA(returned.SpringConstant, 2.0, 1e-9);
        TS_ASSERT_DELTA(returned.DampingConstant, 54.0, 1e-9);
    }

    /**
     * Verify ApplyNScaledInterParameters overwrites mInterEquilibriumDistance and
     * mInterWellDepth using the inter scaling factor (rho).
     */
    void TestApplyNScaledInterParametersOnSemForce()
    {
        SemForce<3> force;
        force.SetInterScalingFactor(3.0);  // different rho than intra

        force.ApplyNScaledInterParameters(27u, 1.0, 9.0, 1.0, 0.7405);

        const double r_eq = 2.0 * std::pow(0.7405 / 27.0, 1.0/3.0);
        const double u0 = 2.0 * r_eq * r_eq / (8.0 * 9.0);  // rho=3, rho^2=9
        TS_ASSERT_DELTA(force.GetInterEquilibriumDistance(), r_eq, 1e-9);
        TS_ASSERT_DELTA(force.GetInterWellDepth(), u0, 1e-9);
        TS_ASSERT_DELTA(force.GetInterScalingFactor(), 3.0, 1e-9);  // rho unchanged
    }

    /**
     * Verify the same Apply methods on SemLinearForce.
     */
    void TestApplyNScaledParametersOnSemLinearForce()
    {
        SemLinearForce<3> force;
        force.SetIntraScalingFactor(5.0);
        // Default eta0 = 1.0, so the returned damping constant should be N = 27
        const SemNScaledParameters returned = force.ApplyNScaledIntraParameters(27u, 1.0, 9.0, 1.0, 0.7405);

        const double r_eq = 2.0 * std::pow(0.7405 / 27.0, 1.0/3.0);
        const double u0 = 2.0 * r_eq * r_eq / (8.0 * 25.0);
        TS_ASSERT_DELTA(force.GetIntraEquilibriumDistance(), r_eq, 1e-9);
        TS_ASSERT_DELTA(force.GetIntraWellDepth(), u0, 1e-9);
        TS_ASSERT_DELTA(returned.DampingConstant, 27.0, 1e-9);
    }
};

#endif /*TESTSEMPARAMETERSCALING_HPP_*/
