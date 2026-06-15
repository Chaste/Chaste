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

#ifndef SEMFORCE_HPP_
#define SEMFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"
#include "SemBasedCellPopulation.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A unified force law for the Subcellular Element Model implementing the
 * modified Morse potential from Sandersius & Newman (2008) Phys. Biol. 5 015002.
 *
 * The potential has the form:
 *   V(r) = u0 * exp(2*rho*(1 - r^2/r_eq^2)) - 2*u0 * exp(rho*(1 - r^2/r_eq^2))
 *
 * This gives a force on node A from node B:
 *   F_A = (4*rho*u0/r_eq^2) * [exp(rho*(1-s)) - exp(2*rho*(1-s))] * (r_B - r_A)
 * where s = r^2/r_eq^2.
 *
 * Separate parameters are used for intra-cellular interactions (nodes in the
 * same element) and inter-cellular interactions (nodes in different elements).
 *
 * The harmonic spring constant near equilibrium is kappa = 8*rho^2*u0/r_eq^2.
 *
 * N-dependent scaling (Sandersius 2008 Section 2):
 *   r_eq(N) = 2*R_cell*(p_3/N)^(1/3)
 *   kappa   = kappa_0 * N^(-1/3) * (1 - lambda*N^(-1/3))
 *   eta     = eta_0 * N
 * These are left as user-configurable via the parameter setters.
 */
template <unsigned DIM>
class SemForce : public AbstractTwoBodyInteractionForce<DIM>
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
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<DIM>>(*this);
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
     * Calculate the force vector on node A due to its interaction with node B,
     * using the full modified Morse potential.
     *
     * F_A = (4*rho*u0/r_eq^2) * [exp(rho*(1-s)) - exp(2*rho*(1-s))] * (r_B - r_A)
     * where s = |r_B - r_A|^2 / r_eq^2.
     *
     * @param rVectorAtoB the vector from node A to node B (r_B - r_A)
     * @param distanceSq the squared distance |r_B - r_A|^2
     * @param u0 the well depth
     * @param rho the scaling factor
     * @param rEq the equilibrium distance
     *
     * @return the force vector on node A
     */
    virtual c_vector<double, DIM> CalculateForceVector(
        const c_vector<double, DIM>& rVectorAtoB,
        double distanceSq,
        double u0,
        double rho,
        double rEq) const;

public:
    /**
     * Constructor.
     */
    SemForce();

    /**
     * Destructor.
     */
    virtual ~SemForce() = default;

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
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemForce)

#endif /*SEMFORCE_HPP_*/
