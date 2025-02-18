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

#include "ImmersedBoundaryKinematicFeedbackForce.hpp"
#include "ImmersedBoundaryEnumerations.hpp"
#include "UblasCustomFunctions.hpp"

template <unsigned DIM>
ImmersedBoundaryKinematicFeedbackForce<DIM>::ImmersedBoundaryKinematicFeedbackForce()
    : AbstractImmersedBoundaryForce<DIM>(),
      mSpringConst(1e3),
      mPreviousLocations()
{
}

template <unsigned DIM>
void ImmersedBoundaryKinematicFeedbackForce<DIM>::AddImmersedBoundaryForceContribution(
    std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
    ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    // Update memory allocation if number of nodes has changed/on first call
    if (mPreviousLocations.size() != rCellPopulation.GetNumNodes())
    {
        mPreviousLocations.resize(rCellPopulation.GetNumNodes());
        UpdatePreviousLocations(rCellPopulation);
    }

    /*
     * Calculate force only if neither node is in a lamina, and if nodes are in
     * different elements.
     */
    auto condition_satisfied = [&rCellPopulation](const std::pair<Node<DIM>*, Node<DIM>*>& pair) -> bool
    {
        // Laminas do not participate in this force class
        if (pair.first->GetRegion() == LAMINA_REGION || pair.second->GetRegion() == LAMINA_REGION)
        {
            return false;
        }
        // This force only acts on nodes in different elements
        if (*(pair.first->ContainingElementsBegin()) == *(pair.second->ContainingElementsBegin()))
        {
            return false;
        }
        // Return true if the nodes are within threshold distance, else false
        auto vec_a2b = rCellPopulation.rGetMesh().GetVectorFromAtoB(pair.first->rGetLocation(), pair.second->rGetLocation());
        return norm_2(vec_a2b) < rCellPopulation.GetInteractionDistance();
    };

    for (auto&& pair : rNodePairs)
    {
        if (condition_satisfied(pair))
        {
            Node<DIM>* p_node_a = pair.first;
            Node<DIM>* p_node_b = pair.second;

            // Get the current and previous displacement between nodes
            auto previous_disp = rCellPopulation.rGetMesh().GetVectorFromAtoB(mPreviousLocations[p_node_a->GetIndex()],
                                                                              mPreviousLocations[p_node_b->GetIndex()]);
            auto current_disp = rCellPopulation.rGetMesh().GetVectorFromAtoB(p_node_a->rGetLocation(),
                                                                             p_node_b->rGetLocation());

            // Calculate the relative velocity and fill in unit_perp, the direction of the force that will act
            c_vector<double, DIM> unit_perp;
            double relative_vel_comp = CalculateRelativeVelocityComponent(previous_disp, current_disp, unit_perp);

            unsigned a_idx = *(p_node_a->ContainingElementsBegin());
            unsigned b_idx = *(p_node_b->ContainingElementsBegin());

            double node_a_elem_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(a_idx, false);
            double node_b_elem_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(b_idx, false);

            double elem_spacing = 0.5 * (node_a_elem_spacing + node_b_elem_spacing);

            double eff_spring_const = mSpringConst * elem_spacing / rCellPopulation.GetIntrinsicSpacing();

            /*
             * We must scale each applied force by a factor of elem_spacing / local spacing, so that forces
             * balance when spread to the grid later (where the multiplicative factor is the local spacing)
             */
            // df = K dx
            c_vector<double, DIM> force = unit_perp * (relative_vel_comp * eff_spring_const);

            c_vector<double, DIM> force_on_b = force * (elem_spacing / node_a_elem_spacing);
            p_node_b->AddAppliedForceContribution(force_on_b);

            c_vector<double, DIM> force_on_a = force * (-1.0 * elem_spacing / node_b_elem_spacing);
            p_node_a->AddAppliedForceContribution(force_on_a);
        }
    }

    UpdatePreviousLocations(rCellPopulation);

    if (this->mAdditiveNormalNoise)
    {
        this->AddNormalNoiseToNodes(rCellPopulation);
    }
}

template <unsigned DIM>
double ImmersedBoundaryKinematicFeedbackForce<DIM>::CalculateRelativeVelocityComponent(
    const c_vector<double, DIM>& rPreviousDisp,
    const c_vector<double, DIM>& rCurrentDisp,
    c_vector<double, DIM>& rUnitPerp)
{
    // Get a unit vector perpendicular to the line joining the nodes at the previous time step
    rUnitPerp = Create_c_vector(-rPreviousDisp[1], rPreviousDisp[0]);
    rUnitPerp /= norm_2(rUnitPerp);

    // Calculate the relative velocity component in the direction of the perpendicular
    return inner_prod(rCurrentDisp / SimulationTime::Instance()->GetTimeStep(), rUnitPerp);
}

template<unsigned DIM>
void ImmersedBoundaryKinematicFeedbackForce<DIM>::UpdatePreviousLocations(
    ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    /*
     * Populate the mPreviousLocations vector with the current location of
     * nodes, so it's ready for next time step.
     */
    for (const auto& p_node : rCellPopulation.rGetMesh().rGetNodes())
    {
        if (p_node->GetRegion() != LAMINA_REGION)
        {
            if (p_node->GetIndex() >= mPreviousLocations.size())
            {
                mPreviousLocations.resize(p_node->GetIndex() + 1);
            }

            mPreviousLocations[p_node->GetIndex()] = p_node->rGetLocation();
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryKinematicFeedbackForce<DIM>::OutputImmersedBoundaryForceParameters(
    out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SpringConstant>" << mSpringConst << "</SpringConstant>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template<unsigned DIM>
double ImmersedBoundaryKinematicFeedbackForce<DIM>::GetSpringConst() const
{
    return mSpringConst;
}

template<unsigned DIM>
void ImmersedBoundaryKinematicFeedbackForce<DIM>::SetSpringConst(
    double springConst)
{
    mSpringConst = springConst;
}

// Explicit instantiation
template class ImmersedBoundaryKinematicFeedbackForce<1>;
template class ImmersedBoundaryKinematicFeedbackForce<2>;
template class ImmersedBoundaryKinematicFeedbackForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryKinematicFeedbackForce)
