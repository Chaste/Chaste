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

#ifndef SEMREGIONALFORCE_HPP_
#define SEMREGIONALFORCE_HPP_

#include <vector>

#include "AbstractTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A linear (harmonic) spring force for Subcellular Element (SEM) simulations in which the spring
 * constant and rest length depend on the region label of the interacting nodes.
 *
 * Mechanically this is identical to SemLinearForce: the force on node A from node B is the
 * harmonic spring
 *   F_A = kappa * (1 - L / |r|) * (r_B - r_A),
 * zero beyond a cut-off distance and zero for coincident nodes. The only difference is that the
 * spring constant kappa and rest length L are chosen per node region rather than from an
 * intra-/inter-cellular parameter set.
 *
 * Each subcellular node carries an unsigned region label (Node::GetRegion()) set by the mesh
 * generator; by default the SemNodeRegion enum labels the interior (0) and boundary/cortex (1),
 * but any number of regions is supported. There is one spring constant and one rest length per
 * region (mSpringConstants and mRestLengths, indexed by region label). For a pair of interacting
 * nodes the effective spring constant and rest length are the arithmetic mean of each node's
 * per-region value, allowing smooth transitions where regions meet.
 */
template<unsigned DIM>
class SemRegionalForce : public AbstractTwoBodyInteractionForce<DIM>
{
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
        archive & mSpringConstants;
        archive & mRestLengths;
        archive & mCutOffDistance;
    }

protected:

    /** Spring constant per node region, indexed by region label. */
    std::vector<double> mSpringConstants;

    /** Rest length per node region, indexed by region label. */
    std::vector<double> mRestLengths;

    /** Cut-off distance beyond which the interaction is zero. */
    double mCutOffDistance;

    /**
     * Calculate the linear (harmonic) spring force vector for a SEM interaction.
     *
     * F_A = springConstant * (1 - restLength / |r|) * (r_B - r_A)
     *
     * Coincident nodes (|r| ~ 0) return the zero vector.
     *
     * @param rVectorAtoB the vector from node A to node B (r_B - r_A)
     * @param distanceSq the squared distance |r_B - r_A|^2
     * @param springConstant the (combined) spring constant
     * @param restLength the (combined) rest length
     *
     * @return the force vector on node A
     */
    c_vector<double, DIM> CalculateForceVector(
        const c_vector<double, DIM>& rVectorAtoB,
        double distanceSq,
        double springConstant,
        double restLength) const;

public:

    /**
     * Constructor. Defaults to two regions (interior, boundary) matching the SemNodeRegion enum.
     */
    SemRegionalForce();

    /**
     * Destructor.
     */
    virtual ~SemRegionalForce() = default;

    /**
     * Overridden AddForceContribution() method.
     *
     * This force is only valid for SEM populations, and requires the per-region spring-constant
     * and rest-length arrays to be non-empty and of equal length. Once these are verified, force
     * application is delegated to the generic two-body interaction machinery in the parent class.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation) override;

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the harmonic spring force between two nodes using region-dependent spring
     * constants and rest lengths (the arithmetic mean of the two nodes' per-region values).
     *
     * Note that this assumes the nodes are neighbours and is called by AddForceContribution().
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
     * Set the per-region spring constants. Its length (which must match the rest-length array's
     * length, checked at force-application time) determines the number of regions.
     *
     * @param rSpringConstants one spring constant per region
     */
    void SetSpringConstants(const std::vector<double>& rSpringConstants);

    /**
     * @return the per-region spring constants
     */
    std::vector<double> GetSpringConstants() const;

    /**
     * Set the per-region rest lengths. Its length (which must match the spring-constant array's
     * length, checked at force-application time) determines the number of regions.
     *
     * @param rRestLengths one rest length per region
     */
    void SetRestLengths(const std::vector<double>& rRestLengths);

    /**
     * @return the per-region rest lengths
     */
    std::vector<double> GetRestLengths() const;

    /**
     * @param cutOffDistance the new cut-off distance
     */
    void SetCutOffDistance(double cutOffDistance);

    /**
     * @return the cut-off distance
     */
    double GetCutOffDistance() const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemRegionalForce)

#endif /*SEMREGIONALFORCE_HPP_*/
