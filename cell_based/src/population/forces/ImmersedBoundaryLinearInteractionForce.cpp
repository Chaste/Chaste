/*

Copyright (c) 2005-2025, University of Oxford.
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

#include "ImmersedBoundaryLinearInteractionForce.hpp"
#include "ImmersedBoundaryEnumerations.hpp"

template <unsigned DIM>
ImmersedBoundaryLinearInteractionForce<DIM>::ImmersedBoundaryLinearInteractionForce()
    : AbstractImmersedBoundaryForce<DIM>(),
       mSpringConst(1e3),
       mRestLength(0.25),
       mLaminaSpringConstMult(1.0),
       mLaminaRestLengthMult(1.0)
{
}

template<unsigned DIM>
ImmersedBoundaryLinearInteractionForce<DIM>::~ImmersedBoundaryLinearInteractionForce()
{
}

template <unsigned DIM>
void ImmersedBoundaryLinearInteractionForce<DIM>::AddImmersedBoundaryForceContribution(
    std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
    ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    for (unsigned pair = 0; pair < rNodePairs.size(); ++pair)
    {
        /*
         * Interactions only exist between pairs of nodes that are not in the
         * same boundary / lamina.
         */
        if (rCellPopulation.rGetMesh().NodesInDifferentElementOrLamina(rNodePairs[pair].first, rNodePairs[pair].second))
        {
            Node<DIM>* p_node_a = rNodePairs[pair].first;
            Node<DIM>* p_node_b = rNodePairs[pair].second;

            c_vector<double, DIM> vec_a2b = rCellPopulation.rGetMesh().GetVectorFromAtoB(p_node_a->rGetLocation(),
                                                                                         p_node_b->rGetLocation());
            double normed_dist = norm_2(vec_a2b);

            // Force non-zero only within interaction distance, by definition
            if (normed_dist < rCellPopulation.GetInteractionDistance())
            {
                // Need the average spacing of the containing element; this depends on whether it's a lamina or element
                bool a_lamina = p_node_a->GetRegion() == LAMINA_REGION;
                bool b_lamina = p_node_b->GetRegion() == LAMINA_REGION;

                unsigned a_idx = *(p_node_a->ContainingElementsBegin());
                unsigned b_idx = *(p_node_b->ContainingElementsBegin());

                double node_a_elem_spacing = a_lamina ? rCellPopulation.rGetMesh().GetAverageNodeSpacingOfLamina(a_idx, false)
                                                      : rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(a_idx, false);

                double node_b_elem_spacing = b_lamina ? rCellPopulation.rGetMesh().GetAverageNodeSpacingOfLamina(b_idx, false)
                                                      : rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(b_idx, false);

                double elem_spacing = 0.5 * (node_a_elem_spacing + node_b_elem_spacing);

                double eff_spring_const = mSpringConst * elem_spacing / rCellPopulation.GetIntrinsicSpacing();
                double eff_rest_length = mRestLength * rCellPopulation.GetInteractionDistance();

                if (a_lamina || b_lamina)
                {
                    eff_spring_const *= mLaminaSpringConstMult;
                    eff_rest_length *= mLaminaRestLengthMult;
                }

                /*
                 * We must scale each applied force by a factor of elem_spacing
                 * / local spacing, so that forces balance when spread to the
                 * grid later (where the multiplicative factor is the local
                 * spacing).
                 */
                vec_a2b *= eff_spring_const * (normed_dist - eff_rest_length) / normed_dist;

                c_vector<double, DIM> force_a2b = vec_a2b * (elem_spacing / node_a_elem_spacing);
                p_node_a->AddAppliedForceContribution(force_a2b);

                c_vector<double, DIM> force_b2a = vec_a2b * (-1.0 * elem_spacing / node_b_elem_spacing);
                p_node_b->AddAppliedForceContribution(force_b2a);
            }
        }
    }

    if (this->mAdditiveNormalNoise)
    {
        this->AddNormalNoiseToNodes(rCellPopulation);
    }
}


template<unsigned DIM>
void ImmersedBoundaryLinearInteractionForce<DIM>::OutputImmersedBoundaryForceParameters(
    out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SpringConstant>" << mSpringConst << "</SpringConstant>\n";
    *rParamsFile << "\t\t\t<RestLength>" << mRestLength << "</RestLength>\n";
    *rParamsFile << "\t\t\t<LaminaSpringConstMult>" << mLaminaSpringConstMult << "</LaminaSpringConstMult>\n";
    *rParamsFile << "\t\t\t<LaminaRestLengthMult>" << mLaminaRestLengthMult << "</LaminaRestLengthMult>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template<unsigned DIM>
double ImmersedBoundaryLinearInteractionForce<DIM>::GetSpringConst() const
{
    return mSpringConst;
}

template<unsigned DIM>
void ImmersedBoundaryLinearInteractionForce<DIM>::SetSpringConst(
    double springConst)
{
    mSpringConst = springConst;
}

template<unsigned DIM>
double ImmersedBoundaryLinearInteractionForce<DIM>::GetRestLength() const
{
    return mRestLength;
}

template<unsigned DIM>
void ImmersedBoundaryLinearInteractionForce<DIM>::SetRestLength(double restLength)
{
    mRestLength = restLength;
}

template<unsigned DIM>
double ImmersedBoundaryLinearInteractionForce<DIM>::GetLaminaSpringConstMult() const
{
    return mLaminaSpringConstMult;
}

template<unsigned DIM>
void ImmersedBoundaryLinearInteractionForce<DIM>::SetLaminaSpringConstMult(
    double laminaSpringConstMult)
{
    mLaminaSpringConstMult = laminaSpringConstMult;
}

template<unsigned DIM>
double ImmersedBoundaryLinearInteractionForce<DIM>::GetLaminaRestLengthMult() const
{
    return mLaminaRestLengthMult;
}

template<unsigned DIM>
void ImmersedBoundaryLinearInteractionForce<DIM>::SetLaminaRestLengthMult(
    double laminaRestLengthMult)
{
    mLaminaRestLengthMult = laminaRestLengthMult;
}

// Explicit instantiation
template class ImmersedBoundaryLinearInteractionForce<1>;
template class ImmersedBoundaryLinearInteractionForce<2>;
template class ImmersedBoundaryLinearInteractionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryLinearInteractionForce)
