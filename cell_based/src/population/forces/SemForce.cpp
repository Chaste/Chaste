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

#include "SemForce.hpp"

template<unsigned DIM>
SemForce<DIM>::SemForce()
   : AbstractTwoBodyInteractionForce<DIM>(),
     mIntraWellDepth(1.0),
     mIntraScalingFactor(5.0),
     mIntraEquilibriumDistance(0.2),
     mIntraCutOffDistance(0.5),
     mInterWellDepth(1.0),
     mInterScalingFactor(5.0),
     mInterEquilibriumDistance(0.3),
     mInterCutOffDistance(0.5)
{
}

template<unsigned DIM>
c_vector<double, DIM> SemForce<DIM>::CalculateForceVector(
    const c_vector<double, DIM>& rVectorAtoB,
    double distanceSq,
    double u0,
    double rho,
    double rEq) const
{
    /*
     * Full modified Morse potential (Sandersius & Newman 2008):
     *   V(r) = u0 * exp(2*rho*(1 - r^2/r_eq^2)) - 2*u0 * exp(rho*(1 - r^2/r_eq^2))
     *
     * Force on A = -grad_A V = (dV/dr) * r_hat, with r_hat = (r_B - r_A)/|r_B - r_A|.
     *
     * Taking the derivative and noting the factor of r in dV/dr cancels with
     * the normalisation of r_hat, we get:
     *   F_A = (4*rho*u0/r_eq^2) * [exp(rho*(1-s)) - exp(2*rho*(1-s))] * (r_B - r_A)
     * where s = r^2/r_eq^2.
     *
     * Sign convention: positive coefficient => force toward B (attractive);
     * negative coefficient => force away from B (repulsive).
     * At s < 1 (r < r_eq): exp(rho*(1-s)) < exp(2*rho*(1-s)), so coeff < 0 => repulsive. Correct.
     * At s > 1 (r > r_eq): exp(rho*(1-s)) > exp(2*rho*(1-s)), so coeff > 0 => attractive. Correct.
     */
    const double rEqSq = rEq * rEq;
    const double s = distanceSq / rEqSq;
    const double exp_rho = std::exp(rho * (1.0 - s));
    const double exp_2rho = exp_rho * exp_rho;  // exp(2*rho*(1-s)) = [exp(rho*(1-s))]^2

    const double coefficient = (4.0 * rho * u0 / rEqSq) * (exp_rho - exp_2rho);

    return coefficient * rVectorAtoB;
}

template<unsigned DIM>
void SemForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if (dynamic_cast<SemBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SemForce is to be used with a SemBasedCellPopulation only");
    }

    AbstractTwoBodyInteractionForce<DIM>::AddForceContribution(rCellPopulation);
}

template<unsigned DIM>
c_vector<double, DIM> SemForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                unsigned nodeBGlobalIndex,
                                                                AbstractCellPopulation<DIM>& rCellPopulation)
{
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    if (dynamic_cast<SemBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SemForce is to be used with a SemBasedCellPopulation only");
    }

    Node<DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    const c_vector<double, DIM>& r_loc_a = p_node_a->rGetLocation();
    const c_vector<double, DIM>& r_loc_b = p_node_b->rGetLocation();
    const c_vector<double, DIM> vec_a_to_b = r_loc_b - r_loc_a;
    const double dist_sq = inner_prod(vec_a_to_b, vec_a_to_b);

    const unsigned elem_a = *p_node_a->rGetContainingElementIndices().begin();
    const unsigned elem_b = *p_node_b->rGetContainingElementIndices().begin();
    const bool same_cell = (elem_a == elem_b);

    const double u0 = same_cell ? mIntraWellDepth : mInterWellDepth;
    const double rho = same_cell ? mIntraScalingFactor : mInterScalingFactor;
    const double rEq = same_cell ? mIntraEquilibriumDistance : mInterEquilibriumDistance;
    const double cutoff = same_cell ? mIntraCutOffDistance : mInterCutOffDistance;

    if (dist_sq >= cutoff * cutoff)
    {
        return zero_vector<double>(DIM);
    }

    return CalculateForceVector(vec_a_to_b, dist_sq, u0, rho, rEq);
}

template<unsigned DIM>
void SemForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<IntraWellDepth>" << mIntraWellDepth << "</IntraWellDepth>\n";
    *rParamsFile << "\t\t\t<IntraScalingFactor>" << mIntraScalingFactor << "</IntraScalingFactor>\n";
    *rParamsFile << "\t\t\t<IntraEquilibriumDistance>" << mIntraEquilibriumDistance << "</IntraEquilibriumDistance>\n";
    *rParamsFile << "\t\t\t<IntraCutOffDistance>" << mIntraCutOffDistance << "</IntraCutOffDistance>\n";
    *rParamsFile << "\t\t\t<InterWellDepth>" << mInterWellDepth << "</InterWellDepth>\n";
    *rParamsFile << "\t\t\t<InterScalingFactor>" << mInterScalingFactor << "</InterScalingFactor>\n";
    *rParamsFile << "\t\t\t<InterEquilibriumDistance>" << mInterEquilibriumDistance << "</InterEquilibriumDistance>\n";
    *rParamsFile << "\t\t\t<InterCutOffDistance>" << mInterCutOffDistance << "</InterCutOffDistance>\n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

// Getters and setters

template<unsigned DIM>
double SemForce<DIM>::GetIntraWellDepth() const { return mIntraWellDepth; }

template<unsigned DIM>
void SemForce<DIM>::SetIntraWellDepth(double wellDepth) { mIntraWellDepth = wellDepth; }

template<unsigned DIM>
double SemForce<DIM>::GetIntraScalingFactor() const { return mIntraScalingFactor; }

template<unsigned DIM>
void SemForce<DIM>::SetIntraScalingFactor(double scalingFactor) { mIntraScalingFactor = scalingFactor; }

template<unsigned DIM>
double SemForce<DIM>::GetIntraEquilibriumDistance() const { return mIntraEquilibriumDistance; }

template<unsigned DIM>
void SemForce<DIM>::SetIntraEquilibriumDistance(double equilibriumDistance) { mIntraEquilibriumDistance = equilibriumDistance; }

template<unsigned DIM>
double SemForce<DIM>::GetIntraCutOffDistance() const { return mIntraCutOffDistance; }

template<unsigned DIM>
void SemForce<DIM>::SetIntraCutOffDistance(double cutOffDistance) { mIntraCutOffDistance = cutOffDistance; }

template<unsigned DIM>
double SemForce<DIM>::GetInterWellDepth() const { return mInterWellDepth; }

template<unsigned DIM>
void SemForce<DIM>::SetInterWellDepth(double wellDepth) { mInterWellDepth = wellDepth; }

template<unsigned DIM>
double SemForce<DIM>::GetInterScalingFactor() const { return mInterScalingFactor; }

template<unsigned DIM>
void SemForce<DIM>::SetInterScalingFactor(double scalingFactor) { mInterScalingFactor = scalingFactor; }

template<unsigned DIM>
double SemForce<DIM>::GetInterEquilibriumDistance() const { return mInterEquilibriumDistance; }

template<unsigned DIM>
void SemForce<DIM>::SetInterEquilibriumDistance(double equilibriumDistance) { mInterEquilibriumDistance = equilibriumDistance; }

template<unsigned DIM>
double SemForce<DIM>::GetInterCutOffDistance() const { return mInterCutOffDistance; }

template<unsigned DIM>
void SemForce<DIM>::SetInterCutOffDistance(double cutOffDistance) { mInterCutOffDistance = cutOffDistance; }

template<unsigned DIM>
SemNScaledParameters SemForce<DIM>::ApplyNScaledIntraParameters(
    unsigned numNodes, double cellRadius, double kappa0, double lambda, double packingDensity, double eta0)
{
    const SemNScaledParameters params = SemComputeNScaledParameters<DIM>(
        numNodes, cellRadius, kappa0, mIntraScalingFactor, lambda, eta0, packingDensity);
    mIntraEquilibriumDistance = params.EquilibriumDistance;
    mIntraWellDepth = params.WellDepth;
    return params;
}

template<unsigned DIM>
SemNScaledParameters SemForce<DIM>::ApplyNScaledInterParameters(
    unsigned numNodes, double cellRadius, double kappa0, double lambda, double packingDensity, double eta0)
{
    const SemNScaledParameters params = SemComputeNScaledParameters<DIM>(
        numNodes, cellRadius, kappa0, mInterScalingFactor, lambda, eta0, packingDensity);
    mInterEquilibriumDistance = params.EquilibriumDistance;
    mInterWellDepth = params.WellDepth;
    return params;
}

// Explicit instantiation
template class SemForce<1>;
template class SemForce<2>;
template class SemForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemForce)
