/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "ImmersedBoundaryKinematicFeedbackForce.hpp"
#include "ImmersedBoundaryEnumerations.hpp"

template <unsigned DIM>
ImmersedBoundaryKinematicFeedbackForce<DIM>::ImmersedBoundaryKinematicFeedbackForce()
        : AbstractImmersedBoundaryForce<DIM>(),
          mSpringConst(1e3),
          mRestLength(0.25),
          mPreviousLocations({})
{
}

template <unsigned DIM>
void ImmersedBoundaryKinematicFeedbackForce<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                                                                       ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    if (mPreviousLocations.empty())
    {
        mPreviousLocations.resize(rCellPopulation.GetNumNodes());
        UpdatePreviousLocations(rCellPopulation);
    }

    double dt = SimulationTime::Instance()->GetTimeStep();

    // Calculate force only if neither node is in a lamina, and if nodes are in different elements
    auto condition_satisfied = [](const std::pair<Node<DIM>*, Node<DIM>*>& pair) -> bool
    {
        if (pair.first->GetRegion() == LAMINA_REGION || pair.second->GetRegion() == LAMINA_REGION)
        {
            return false;
        }
        return *(pair.first->ContainingElementsBegin()) != *(pair.second->ContainingElementsBegin());
    };

    for (auto && pair : rNodePairs)
    {
        if (condition_satisfied(pair))
        {
            Node<DIM>* p_node_a = pair.first;
            Node<DIM>* p_node_b = pair.second;

            // Get unit vector between nodes
            auto vec_a2b = rCellPopulation.rGetMesh().GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());
            vec_a2b /= norm_2(vec_a2b);

            auto disp_a = p_node_a->rGetLocation() - mPreviousLocations[p_node_a->GetIndex()];
            auto disp_b = p_node_b->rGetLocation() - mPreviousLocations[p_node_b->GetIndex()];

            auto velocity = (disp_a - disp_b) / dt;
        }
    }

    for (unsigned pair = 0; pair < rNodePairs.size(); pair++)
    {
        // Interactions only exist between pairs of nodes that are not in the same boundary / lamina
        if (true)
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


                /*
                 * We must scale each applied force by a factor of elem_spacing / local spacing, so that forces
                 * balance when spread to the grid later (where the multiplicative factor is the local spacing)
                 */
                vec_a2b *= eff_spring_const * (normed_dist - eff_rest_length) / normed_dist;

                c_vector<double, DIM> force_a2b = vec_a2b * (elem_spacing / node_a_elem_spacing);
                p_node_a->AddAppliedForceContribution(force_a2b);

                c_vector<double, DIM> force_b2a = vec_a2b * (-1.0 * elem_spacing / node_b_elem_spacing);
                p_node_b->AddAppliedForceContribution(force_b2a);
            }
        }
    }

    UpdatePreviousLocations(rCellPopulation);

    if (this->mAdditiveNormalNoise)
    {
        this->AddNormalNoiseToNodes(rCellPopulation);
    }
}

template<unsigned DIM>
void ImmersedBoundaryKinematicFeedbackForce<DIM>::UpdatePreviousLocations(ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    for (auto && p_node : rCellPopulation.rGetMesh().rGetNodes())
    {
        if (p_node->GetRegion() != LAMINA_REGION)
        {
            if (p_node->GetIndex() > mPreviousLocations.size())
            {
                mPreviousLocations.resize(p_node->GetIndex());
            }

            mPreviousLocations[p_node->GetIndex()] = p_node->rGetLocation();
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryKinematicFeedbackForce<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SpringConstant>" << mSpringConst << "</SpringConstant>\n";
    *rParamsFile << "\t\t\t<RestLength>" << mRestLength << "</RestLength>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template<unsigned DIM>
double ImmersedBoundaryKinematicFeedbackForce<DIM>::GetSpringConst() const
{
    return mSpringConst;
}

template<unsigned DIM>
void ImmersedBoundaryKinematicFeedbackForce<DIM>::SetSpringConst(double springConst)
{
    mSpringConst = springConst;
}

template<unsigned DIM>
double ImmersedBoundaryKinematicFeedbackForce<DIM>::GetRestLength() const
{
    return mRestLength;
}

template<unsigned DIM>
void ImmersedBoundaryKinematicFeedbackForce<DIM>::SetRestLength(double restLength)
{
    mRestLength = restLength;
}

// Explicit instantiation
template class ImmersedBoundaryKinematicFeedbackForce<1>;
template class ImmersedBoundaryKinematicFeedbackForce<2>;
template class ImmersedBoundaryKinematicFeedbackForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryKinematicFeedbackForce)
