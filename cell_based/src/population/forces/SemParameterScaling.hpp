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

#ifndef SEMPARAMETERSCALING_HPP_
#define SEMPARAMETERSCALING_HPP_

#include <cmath>
#include <cassert>

/**
 * Parameters computed by N-dependent scaling for SEM simulations.
 *
 * Declared as a class with accessors rather than an aggregate of public members
 * because cppwg cannot bind public data members, so a plain struct would reach
 * Python with no readable fields.
 *
 * @see SemComputeNScaledParameters
 */
class SemNScaledParameters
{
private:

    /** Equilibrium distance r_eq = 2 R_cell (p/N)^{1/DIM}. */
    double mEquilibriumDistance;

    /** Morse well depth u₀ = κ r_eq² / (8 ρ²), derived from the spring constant and rho. */
    double mWellDepth;

    /**
     * Harmonic spring constant κ = κ₀ N^{-1/DIM} (1 − λ N^{-1/DIM}).
     * Directly usable as the spring constant for SemLinearForce.
     */
    double mSpringConstant;

    /**
     * Damping constant η = η₀ N.
     * Pass to SemBasedCellPopulation::SetDampingConstantNormal().
     */
    double mDampingConstant;

public:

    /**
     * Constructor.
     *
     * @param equilibriumDistance  equilibrium distance r_eq
     * @param wellDepth  Morse well depth u₀
     * @param springConstant  harmonic spring constant κ
     * @param dampingConstant  damping constant η
     */
    SemNScaledParameters(double equilibriumDistance,
                         double wellDepth,
                         double springConstant,
                         double dampingConstant)
        : mEquilibriumDistance(equilibriumDistance),
          mWellDepth(wellDepth),
          mSpringConstant(springConstant),
          mDampingConstant(dampingConstant)
    {
    }

    /** @return the equilibrium distance r_eq. */
    double GetEquilibriumDistance() const
    {
        return mEquilibriumDistance;
    }

    /** @return the Morse well depth u₀. */
    double GetWellDepth() const
    {
        return mWellDepth;
    }

    /** @return the harmonic spring constant κ. */
    double GetSpringConstant() const
    {
        return mSpringConstant;
    }

    /** @return the damping constant η. */
    double GetDampingConstant() const
    {
        return mDampingConstant;
    }
};

/**
 * Compute N-dependent Morse/spring parameters for SEM simulations.
 *
 * Implements the scaling laws from Sandersius & Newman (2008) Phys. Biol. 5 015002,
 * Section 2:
 *
 *   r_eq(N)  = 2 R_cell (p/N)^{1/DIM}
 *   κ(N)     = κ₀ N^{−1/DIM} (1 − λ N^{−1/DIM})
 *   u₀       = κ r_eq² / (8 ρ²)           [from κ = 8ρ²u₀/r_eq²]
 *   η        = η₀ N
 *
 * The exponent 1/DIM generalises the 3D paper formulation to 2D and 1D.
 *
 * Default packing densities (can be overridden via packingDensity):
 *   DIM=1:  1.0        (linear)
 *   DIM=2:  π/(2√3) ≈ 0.9069  (hexagonal close packing)
 *   DIM=3:  π/(3√2) ≈ 0.7405  (FCC close packing)
 *
 * The Sandersius paper uses p₃ ≈ 0.64 for random close packing in 3D.
 * Pass packingDensity=0.64 explicitly to reproduce those simulations.
 *
 * The cut-off distance is not set by this function; a typical choice is
 * 1.5–2.5 × EquilibriumDistance, depending on the simulation.
 *
 * @tparam DIM  spatial dimension (sets the scaling exponent to 1/DIM)
 * @param numNodes         number of subcellular nodes N
 * @param cellRadius       cell radius R_cell (same length units as the force's r_eq)
 * @param kappa0           reference spring constant κ₀
 * @param scalingFactor    Morse steepness ρ (used to derive u₀; must equal the force's rho)
 * @param lambda           correction factor λ (default 0; set >0 for softer large-N scaling)
 * @param eta0             reference damping constant per element η₀ (default 1)
 * @param packingDensity   element packing density p (0 → dimension-specific default above)
 * @return                 computed parameter set
 */
template<unsigned DIM>
SemNScaledParameters SemComputeNScaledParameters(
    unsigned numNodes,
    double cellRadius,
    double kappa0,
    double scalingFactor,
    double lambda = 0.0,
    double eta0 = 1.0,
    double packingDensity = 0.0)
{
    assert(numNodes >= 1u);
    assert(cellRadius > 0.0);
    assert(kappa0 > 0.0);
    assert(scalingFactor > 0.0);
    assert(eta0 > 0.0);

    if (packingDensity == 0.0)
    {
        if constexpr (DIM == 1)
            packingDensity = 1.0;
        else if constexpr (DIM == 2)
            packingDensity = 0.9069;  // pi / (2*sqrt(3)), hexagonal close packing
        else
            packingDensity = 0.7405;  // pi / (3*sqrt(2)), FCC close packing
    }

    assert(packingDensity > 0.0);

    const double n_d = static_cast<double>(numNodes);
    const double exponent = 1.0 / static_cast<double>(DIM);

    const double r_eq = 2.0 * cellRadius * std::pow(packingDensity / n_d, exponent);

    const double n_factor = std::pow(n_d, -exponent);
    const double kappa = kappa0 * n_factor * (1.0 - lambda * n_factor);

    const double u0 = kappa * r_eq * r_eq / (8.0 * scalingFactor * scalingFactor);

    const double eta = eta0 * n_d;

    return SemNScaledParameters(r_eq, u0, kappa, eta);
}

#endif // SEMPARAMETERSCALING_HPP_
