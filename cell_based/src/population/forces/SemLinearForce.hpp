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

#ifndef SEMLINEARFORCE_HPP_
#define SEMLINEARFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"
#include "SemBasedCellPopulation.hpp"
#include "SemParameterScaling.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A linear (harmonic) approximation to the SEM Morse force from
 * Sandersius & Newman (2008) Phys. Biol. 5 015002.
 *
 * Uses the same parameters as SemForce (u0, rho, r_eq) but replaces
 * the full Morse potential with its harmonic approximation near
 * equilibrium.
 *
 * The spring constant is derived from the Morse parameters:
 *   kappa = 8 * rho^2 * u0 / r_eq^2
 *
 * The resulting linear force on node A from node B is:
 *   F_A = kappa * (1 - r_eq / |r|) * (r_B - r_A)
 *
 * This should agree well with the full Morse force for small
 * deviations about the equilibrium distance r_eq.
 */
template<unsigned DIM>
class SemLinearForce : public AbstractTwoBodyInteractionForce<DIM>
{
    friend class TestForces;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<DIM> >(*this);
        archive & mIntraWellDepth;
        archive & mIntraScalingFactor;
        archive & mIntraEquilibriumDistance;
        archive & mIntraCutOffDistance;
        archive & mInterWellDepth;
        archive & mInterScalingFactor;
        archive & mInterEquilibriumDistance;
        archive & mInterCutOffDistance;
    }

protected:
    /** Depth of the intra-cellular Morse potential well (u0). */
    double mIntraWellDepth;

    /** Scaling factor for the intra-cellular potential (rho). */
    double mIntraScalingFactor;

    /** Equilibrium distance for the intra-cellular potential (r_eq). */
    double mIntraEquilibriumDistance;

    /** Cut-off distance beyond which intra-cellular interactions are zero. */
    double mIntraCutOffDistance;

    /** Depth of the inter-cellular Morse potential well (u0). */
    double mInterWellDepth;

    /** Scaling factor for the inter-cellular potential (rho). */
    double mInterScalingFactor;

    /** Equilibrium distance for the inter-cellular potential (r_eq). */
    double mInterEquilibriumDistance;

    /** Cut-off distance beyond which inter-cellular interactions are zero. */
    double mInterCutOffDistance;

    /**
     * Calculate the linear (harmonic) spring force vector for a SEM interaction.
     *
     * kappa = 8 * rho^2 * u0 / r_eq^2
     * F_A = kappa * (1 - r_eq / |r|) * (r_B - r_A)
     *
     * @param rVectorAtoB the vector from node A to node B (r_B - r_A)
     * @param distanceSq the squared distance |r_B - r_A|^2
     * @param u0 the well depth
     * @param rho the scaling factor
     * @param rEq the equilibrium distance
     *
     * @return the force vector on node A
     */
    c_vector<double, DIM> CalculateForceVector(
        const c_vector<double, DIM>& rVectorAtoB,
        double distanceSq,
        double u0,
        double rho,
        double rEq) const;

public:

    /**
     * Constructor.
     */
    SemLinearForce();

    /**
     * Destructor.
     */
    virtual ~SemLinearForce() = default;

    /**
     * Overridden AddForceContribution() method.
     *
     * This force is only valid for SEM populations. Once the population type
     * has been verified, force application is delegated to the generic
     * two-body interaction machinery in the parent class.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation) override;

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the linearised SEM interaction between two nodes, using the
     * intra- or inter-cellular parameter set according to element membership.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                     unsigned nodeBGlobalIndex,
                                                     AbstractCellPopulation<DIM>& rCellPopulation) override;

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile) override;

    /**
     * @return the intra-cellular well depth
     */
    double GetIntraWellDepth() const;

    /**
     * @param wellDepth the new intra-cellular well depth
     */
    void SetIntraWellDepth(double wellDepth);

    /**
     * @return the intra-cellular scaling factor
     */
    double GetIntraScalingFactor() const;

    /**
     * @param scalingFactor the new intra-cellular scaling factor
     */
    void SetIntraScalingFactor(double scalingFactor);

    /**
     * @return the intra-cellular equilibrium distance
     */
    double GetIntraEquilibriumDistance() const;

    /**
     * @param equilibriumDistance the new intra-cellular equilibrium distance
     */
    void SetIntraEquilibriumDistance(double equilibriumDistance);

    /**
     * @return the intra-cellular cut-off distance
     */
    double GetIntraCutOffDistance() const;

    /**
     * @param cutOffDistance the new intra-cellular cut-off distance
     */
    void SetIntraCutOffDistance(double cutOffDistance);

    /**
     * @return the inter-cellular well depth
     */
    double GetInterWellDepth() const;

    /**
     * @param wellDepth the new inter-cellular well depth
     */
    void SetInterWellDepth(double wellDepth);

    /**
     * @return the inter-cellular scaling factor
     */
    double GetInterScalingFactor() const;

    /**
     * @param scalingFactor the new inter-cellular scaling factor
     */
    void SetInterScalingFactor(double scalingFactor);

    /**
     * @return the inter-cellular equilibrium distance
     */
    double GetInterEquilibriumDistance() const;

    /**
     * @param equilibriumDistance the new inter-cellular equilibrium distance
     */
    void SetInterEquilibriumDistance(double equilibriumDistance);

    /**
     * @return the inter-cellular cut-off distance
     */
    double GetInterCutOffDistance() const;

    /**
     * @param cutOffDistance the new inter-cellular cut-off distance
     */
    void SetInterCutOffDistance(double cutOffDistance);

    /**
     * Apply N-dependent scaling (Sandersius 2008 Section 2) to the intra-cellular
     * parameters, using the current mIntraScalingFactor (rho). Sets
     * mIntraEquilibriumDistance and mIntraWellDepth.
     *
     * The full scaled parameter set is returned so the caller can apply the
     * consistent damping constant eta = eta0*N to the cell population, e.g.
     * rPopulation.SetDampingConstantNormal(ApplyNScaledIntraParameters(...).DampingConstant).
     * The damping constant is not applied automatically, but returning it here
     * avoids recomputing (and possibly mismatching) the scaling.
     *
     * @param numNodes       number of subcellular nodes N
     * @param cellRadius     cell radius R_cell
     * @param kappa0         reference spring constant κ₀
     * @param lambda         correction factor λ (default 0)
     * @param packingDensity packing density p (0 → dimension-specific default)
     * @param eta0           reference damping constant per element η₀ (default 1)
     *
     * @return the full N-scaled parameter set (including the damping constant)
     */
    SemNScaledParameters ApplyNScaledIntraParameters(unsigned numNodes,
                                                     double cellRadius,
                                                     double kappa0,
                                                     double lambda = 0.0,
                                                     double packingDensity = 0.0,
                                                     double eta0 = 1.0);

    /**
     * Apply N-dependent scaling (Sandersius 2008 Section 2) to the inter-cellular
     * parameters, using the current mInterScalingFactor (rho). Sets
     * mInterEquilibriumDistance and mInterWellDepth.
     *
     * As with ApplyNScaledIntraParameters(), the full scaled parameter set is
     * returned so the caller can apply the consistent damping constant to the
     * cell population.
     *
     * @param numNodes       number of subcellular nodes N
     * @param cellRadius     cell radius R_cell
     * @param kappa0         reference spring constant κ₀
     * @param lambda         correction factor λ (default 0)
     * @param packingDensity packing density p (0 → dimension-specific default)
     * @param eta0           reference damping constant per element η₀ (default 1)
     *
     * @return the full N-scaled parameter set (including the damping constant)
     */
    SemNScaledParameters ApplyNScaledInterParameters(unsigned numNodes,
                                                     double cellRadius,
                                                     double kappa0,
                                                     double lambda = 0.0,
                                                     double packingDensity = 0.0,
                                                     double eta0 = 1.0);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemLinearForce)

#endif /*SEMLINEARFORCE_HPP_*/
