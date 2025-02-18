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

#include "ImmersedBoundaryLinearDifferentialAdhesionForce.hpp"

#include "CellLabel.hpp"

template <unsigned DIM>
ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::ImmersedBoundaryLinearDifferentialAdhesionForce()
    : AbstractImmersedBoundaryForce<DIM>(),
      mLabelledCellToLabelledCellSpringConst(1e3),
      mLabelledCellToCellSpringConst(1e3),
      mCellToCellSpringConst(1e3),
      mRestLength(0.25)
{
}

template <unsigned DIM>
void ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::AddImmersedBoundaryForceContribution(
    std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
    ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    if (rCellPopulation.rGetMesh().GetNumLaminas() > 0u)
    {
        EXCEPTION("This force class cannot be used in the presence of lamina elements.");
    }

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
            const double normed_dist = norm_2(vec_a2b);

            // Force non-zero only within interaction distance, by definition
            if (normed_dist < rCellPopulation.GetInteractionDistance())
            {
                const unsigned a_elem_idx = *(p_node_a->ContainingElementsBegin());
                const unsigned b_elem_idx = *(p_node_b->ContainingElementsBegin());

                bool a_cell_labelled = rCellPopulation.GetCellUsingLocationIndex(a_elem_idx)->template HasCellProperty<CellLabel>();
                bool b_cell_labelled = rCellPopulation.GetCellUsingLocationIndex(b_elem_idx)->template HasCellProperty<CellLabel>();

                double node_a_elem_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(a_elem_idx, false);
                double node_b_elem_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(b_elem_idx, false);

                double elem_spacing = 0.5 * (node_a_elem_spacing + node_b_elem_spacing);

                // Here we take account of whether the cells are labelled
                double eff_spring_const = elem_spacing / rCellPopulation.GetIntrinsicSpacing();
                if (a_cell_labelled && b_cell_labelled)
                {
                    eff_spring_const *= mLabelledCellToLabelledCellSpringConst;
                }
                else if (a_cell_labelled || b_cell_labelled)
                {
                    eff_spring_const *= mLabelledCellToCellSpringConst;
                }
                else
                {
                    eff_spring_const *= mCellToCellSpringConst;
                }

                const double eff_rest_length = mRestLength * rCellPopulation.GetInteractionDistance();

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
}

template<unsigned DIM>
void ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::OutputImmersedBoundaryForceParameters(
    out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<LabelledCellToLabelledCellSpringConst>" << mLabelledCellToLabelledCellSpringConst << "</LabelledCellToLabelledCellSpringConst>\n";
    *rParamsFile << "\t\t\t<LabelledCellToCellSpringConst>" << mLabelledCellToCellSpringConst << "</LabelledCellToCellSpringConst>\n";
    *rParamsFile << "\t\t\t<CellToCellSpringConst>" << mCellToCellSpringConst << "</CellToCellSpringConst>\n";
    *rParamsFile << "\t\t\t<RestLength>" << mRestLength << "</RestLength>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template<unsigned DIM>
double ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::GetRestLength() const
{
    return mRestLength;
}

template<unsigned DIM>
void ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::SetRestLength(
    double restLength)
{
    mRestLength = restLength;
}

template<unsigned int DIM>
double ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::GetLabelledCellToLabelledCellSpringConst() const
{
    return mLabelledCellToLabelledCellSpringConst;
}

template<unsigned int DIM>
void ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::SetLabelledCellToLabelledCellSpringConst(
        double labelledCellToLabelledCellSpringConst)
{
    mLabelledCellToLabelledCellSpringConst = labelledCellToLabelledCellSpringConst;
}

template<unsigned int DIM>
double ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::GetLabelledCellToCellSpringConst() const
{
    return mLabelledCellToCellSpringConst;
}

template<unsigned int DIM>
void ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::SetLabelledCellToCellSpringConst(
        double labelledCellToCellSpringConst)
{
    mLabelledCellToCellSpringConst = labelledCellToCellSpringConst;
}

template<unsigned int DIM>
double ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::GetCellToCellSpringConst() const
{
    return mCellToCellSpringConst;
}

template<unsigned int DIM>
void ImmersedBoundaryLinearDifferentialAdhesionForce<DIM>::SetCellToCellSpringConst(
    double cellToCellSpringConst)
{
    mCellToCellSpringConst = cellToCellSpringConst;
}

// Explicit instantiation
template class ImmersedBoundaryLinearDifferentialAdhesionForce<1>;
template class ImmersedBoundaryLinearDifferentialAdhesionForce<2>;
template class ImmersedBoundaryLinearDifferentialAdhesionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryLinearDifferentialAdhesionForce)
