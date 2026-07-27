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

#include "SemRegionalForce.hpp"

#include "SemBasedCellPopulation.hpp"

#include <cmath>

template<unsigned DIM>
SemRegionalForce<DIM>::SemRegionalForce()
   : AbstractTwoBodyInteractionForce<DIM>(),
     mSpringConstants{1.0, 2.0},
     mRestLengths{0.2, 0.15},
     mCutOffDistance(0.5)
{
}

template<unsigned DIM>
c_vector<double, DIM> SemRegionalForce<DIM>::CalculateForceVector(const c_vector<double, DIM>& rVectorAtoB,
                                                                  double distanceSq,
                                                                  double springConstant,
                                                                  double restLength) const
{
    /*
     * Harmonic spring force (identical in form to SemLinearForce):
     *   F_A = springConstant * (r - restLength) * r_hat
     *       = springConstant * (1 - restLength / r) * (r_B - r_A)
     * where r = |r_B - r_A| and r_hat = (r_B - r_A) / r.
     *
     *   r < restLength => repulsive (force away from B)
     *   r > restLength => attractive (force toward B)
     */
    const double distance = std::sqrt(distanceSq);

    // Guard against division by zero if nodes are coincident
    if (distance < 1e-15)
    {
        return zero_vector<double>(DIM);
    }

    const double coefficient = springConstant * (1.0 - restLength / distance);

    return coefficient * rVectorAtoB;
}

template<unsigned DIM>
void SemRegionalForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if (dynamic_cast<SemBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SemRegionalForce is to be used with a SemBasedCellPopulation only");
    }

    if (mSpringConstants.empty() || mSpringConstants.size() != mRestLengths.size())
    {
        EXCEPTION("SemRegionalForce: the spring-constant and rest-length arrays must be non-empty "
                  "and of equal length (one entry per region)");
    }

    AbstractTwoBodyInteractionForce<DIM>::AddForceContribution(rCellPopulation);
}

template<unsigned DIM>
c_vector<double, DIM> SemRegionalForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                        unsigned nodeBGlobalIndex,
                                                                        AbstractCellPopulation<DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    const c_vector<double, DIM>& r_loc_a = p_node_a->rGetLocation();
    const c_vector<double, DIM>& r_loc_b = p_node_b->rGetLocation();
    const c_vector<double, DIM> vec_a_to_b = r_loc_b - r_loc_a;
    const double dist_sq = inner_prod(vec_a_to_b, vec_a_to_b);

    if (dist_sq >= mCutOffDistance * mCutOffDistance)
    {
        return zero_vector<double>(DIM);
    }

    // Every node region must have a configured spring constant and rest length (the array sizes
    // are validated once per step in AddForceContribution()).
    const unsigned region_a = p_node_a->GetRegion();
    const unsigned region_b = p_node_b->GetRegion();
    assert(region_a < mSpringConstants.size());
    assert(region_b < mSpringConstants.size());

    const double spring_const = 0.5 * (mSpringConstants[region_a] + mSpringConstants[region_b]);
    const double rest_length = 0.5 * (mRestLengths[region_a] + mRestLengths[region_b]);

    return CalculateForceVector(vec_a_to_b, dist_sq, spring_const, rest_length);
}

template<unsigned DIM>
void SemRegionalForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    for (unsigned region = 0; region < mSpringConstants.size(); ++region)
    {
        *rParamsFile << "\t\t\t<SpringConstant_" << region << ">" << mSpringConstants[region] << "</SpringConstant_" << region << ">\n";
    }
    for (unsigned region = 0; region < mRestLengths.size(); ++region)
    {
        *rParamsFile << "\t\t\t<RestLength_" << region << ">" << mRestLengths[region] << "</RestLength_" << region << ">\n";
    }
    *rParamsFile << "\t\t\t<CutOffDistance>" << mCutOffDistance << "</CutOffDistance>\n";

    AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

template<unsigned DIM>
void SemRegionalForce<DIM>::SetSpringConstants(const std::vector<double>& rSpringConstants)
{
    mSpringConstants = rSpringConstants;
}

template<unsigned DIM>
std::vector<double> SemRegionalForce<DIM>::GetSpringConstants() const
{
    return mSpringConstants;
}

template<unsigned DIM>
void SemRegionalForce<DIM>::SetRestLengths(const std::vector<double>& rRestLengths)
{
    mRestLengths = rRestLengths;
}

template<unsigned DIM>
std::vector<double> SemRegionalForce<DIM>::GetRestLengths() const
{
    return mRestLengths;
}

template<unsigned DIM>
void SemRegionalForce<DIM>::SetCutOffDistance(double cutOffDistance)
{
    mCutOffDistance = cutOffDistance;
}

template<unsigned DIM>
double SemRegionalForce<DIM>::GetCutOffDistance() const
{
    return mCutOffDistance;
}

// Explicit instantiation
template class SemRegionalForce<1>;
template class SemRegionalForce<2>;
template class SemRegionalForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemRegionalForce)
