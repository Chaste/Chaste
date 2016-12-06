/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "ImmersedBoundaryMembraneElasticityForce.hpp"

template<unsigned DIM>
ImmersedBoundaryMembraneElasticityForce<DIM>::ImmersedBoundaryMembraneElasticityForce()
    : AbstractImmersedBoundaryForce<DIM>(),
      mpMesh(NULL),
      mSpringConstant(1e6),
      mRestLengthMultiplier(0.5),
      mBasementSpringConstantModifier(1.0),
      mBasementRestLengthModifier(1.0)
{
}

template<unsigned DIM>
ImmersedBoundaryMembraneElasticityForce<DIM>::~ImmersedBoundaryMembraneElasticityForce()
{
}

template<unsigned DIM>
double ImmersedBoundaryMembraneElasticityForce<DIM>::GetApicalLengthForElement(unsigned elemIndex)
{
    assert(mpMesh != NULL);

    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mReferenceLocationInAttributesVector;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location];
}

template<unsigned DIM>
double ImmersedBoundaryMembraneElasticityForce<DIM>::GetBasalLengthForElement(unsigned elemIndex)
{
    assert(mpMesh != NULL);

    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mReferenceLocationInAttributesVector + 1;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location];
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
        ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    // This is only run once, on the first pass, and performs a set-up of all necessary class members
    if (mpMesh == NULL)
    {
        mpMesh = &(rCellPopulation.rGetMesh());

        // Verify whether each element has four corners tagged
        unsigned num_corners = mpMesh->GetElement(0)->rGetCornerNodes().size();
        for (unsigned elem_idx = 1; elem_idx < mpMesh->GetNumElements(); elem_idx++)
        {
            if (num_corners != mpMesh->GetElement(elem_idx)->rGetCornerNodes().size())
            {
                EXCEPTION("All elements must have the same number of corners to use this force class.");
            }
        }

        mElementsHaveCorners = (num_corners == 4);

        // If each element has four corners tagged, we set up apical and basal lengths
        if (mElementsHaveCorners)
        {
            // First verify that all elements have the same number of attributes
            mReferenceLocationInAttributesVector = mpMesh->GetElement(0)->GetNumElementAttributes();
            for (unsigned elem_idx = 1; elem_idx < mpMesh->GetNumElements(); elem_idx++)
            {
                if (mReferenceLocationInAttributesVector != mpMesh->GetElement(elem_idx)->GetNumElementAttributes())
                {
                    EXCEPTION("All elements must have the same number of attributes to use this force class.");
                }
            }

            /*
             * We keep track of the initial size of the apical and basal sides.  This will be the initial width of the
             * element.
             *
             * The two element attributes store the apical and basal lengths, respectively.
             *
             * Attribute i:   Apical length
             *           i+1: Basal length
             */
            TagApicalAndBasalLengths();
        }
    }

    /*
     * The following is called each time step, and calculates the forces acting on each element and lamina
     */

    // Used in the calculation of the spring constant
    mIntrinsicSpacingSquared = rCellPopulation.GetIntrinsicSpacing() * rCellPopulation.GetIntrinsicSpacing();

    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = mpMesh->GetElementIteratorBegin();
         elem_it != mpMesh->GetElementIteratorEnd();
         ++elem_it)
    {
        CalculateForcesOnElement(*elem_it);
    }

    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryLaminaIterator lam_it = mpMesh->GetLaminaIteratorBegin();
         lam_it != mpMesh->GetLaminaIteratorEnd();
         ++lam_it)
    {
        CalculateForcesOnElement(*lam_it);
    }
}

template<unsigned DIM>
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::CalculateForcesOnElement(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>& rElement)
{
    // Get index and number of nodes of current element
    unsigned elem_idx = rElement.GetIndex();
    unsigned num_nodes = rElement.GetNumNodes();

    // Helper variables
    double normed_dist;
    c_vector<double, DIM> aggregate_force;

    // Make a vector to store the force on node i+1 from node i
    std::vector<c_vector<double, DIM> > elastic_force_to_next_node(num_nodes);

    /*
     * Get the node spacing ratio for this element.  The rest length and spring constant are derived from this
     * characteristic length.
     *
     * The spring constant is derived with reference to the intrinsic spacing, so that with different node spacings
     * the user-defined parameters do not have to be updated.
     *
     * The correct factor to increase the spring constant by is (intrinsic spacing / node_spacing)^2.  One factor
     * takes into account the energy considerations of the elastic springs, and the other takes account of the
     * factor of node_spacing used in discretising the force relation.
     */

    double node_spacing;
    double spring_constant;
    double rest_length;

    // Determine if we're in a lamina or not
    if(ELEMENT_DIM < SPACE_DIM)
    {
        node_spacing = mpMesh->GetAverageNodeSpacingOfLamina(elem_idx, false);

        spring_constant = mBasementSpringConstantModifier * mSpringConstant * mIntrinsicSpacingSquared / (node_spacing * node_spacing);
        rest_length = mBasementRestLengthModifier * mRestLengthMultiplier * node_spacing;
    }
    else // regular element
    {
        node_spacing = mpMesh->GetAverageNodeSpacingOfElement(elem_idx, false);

        spring_constant = mSpringConstant * mIntrinsicSpacingSquared / (node_spacing * node_spacing);
        rest_length = mRestLengthMultiplier * node_spacing;
    }

    // Loop over nodes and calculate the force exerted on node i+1 by node i
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // Index of the next node, calculated modulo number of nodes in this element
        unsigned next_idx = (node_idx + 1) % num_nodes;

        double modified_spring_constant = spring_constant;
        double modified_rest_length = rest_length;

        // Hooke's law linear spring force
        elastic_force_to_next_node[node_idx] = mpMesh->GetVectorFromAtoB(rElement.GetNodeLocation(node_idx), rElement.GetNodeLocation(next_idx));
        normed_dist = norm_2(elastic_force_to_next_node[node_idx]);
        elastic_force_to_next_node[node_idx] *= modified_spring_constant * (normed_dist - modified_rest_length) / normed_dist;
    }

    // Add the contributions of springs adjacent to each node
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // Get index of previous node
        unsigned prev_idx = (node_idx + num_nodes - 1) % num_nodes;

        aggregate_force = elastic_force_to_next_node[node_idx] - elastic_force_to_next_node[prev_idx];

        // Add the aggregate force contribution to the node
        rElement.GetNode(node_idx)->AddAppliedForceContribution(aggregate_force);
    }

    ///\todo Why is this code commented out?
    // If corners are present, we add on the additional functionality
//        if (mElementsHaveCorners)
//        {
//            // Add force contributions from apical and basal surfaces
//            if (elem_idx != mpMesh->GetMembraneIndex())
//            {
//                std::vector<Node<DIM> *> r_corners = rElement.rGetCornerNodes();
//
//                // Apical surface
//                c_vector<double, DIM> apical_force = mpMesh->GetVectorFromAtoB(r_corners[0]->rGetLocation(),
//                                                                               r_corners[1]->rGetLocation());
//                normed_dist = norm_2(apical_force);
//                apical_force *=
//                        0.0 * mSpringConstant * (normed_dist - GetApicalLengthForElement(elem_idx)) / normed_dist;
//
//                r_corners[0]->AddAppliedForceContribution(apical_force);
//                apical_force *= -1.0;
//                r_corners[1]->AddAppliedForceContribution(apical_force);
//
//                // Basal surface
//                c_vector<double, DIM> basal_force = mpMesh->GetVectorFromAtoB(r_corners[3]->rGetLocation(),
//                                                                              r_corners[2]->rGetLocation());
//
//                normed_dist = norm_2(basal_force);
//                basal_force *=
//                        0.0 * mSpringConstant * (normed_dist - GetBasalLengthForElement(elem_idx)) / normed_dist;
//
//                r_corners[3]->AddAppliedForceContribution(basal_force);
//                basal_force *= -1.0;
//                r_corners[2]->AddAppliedForceContribution(basal_force);
//            }
//        }
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::TagApicalAndBasalLengths()
{
    assert(mpMesh != NULL);

    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = mpMesh->GetElementIteratorBegin();
         elem_it != mpMesh->GetElementIteratorEnd();
         ++elem_it)
    {
        // Elements start roughly rectangular, so the correct apical and basal lengths are the width of the element

        unsigned half_way = unsigned(elem_it->GetNumNodes() / 2.0);

        double elem_width = fabs( mpMesh->GetVectorFromAtoB(elem_it->GetNode(0)->rGetLocation(),
                                                            elem_it->GetNode(half_way)->rGetLocation())[0] );

        elem_it->AddElementAttribute(elem_width);
        elem_it->AddElementAttribute(elem_width);
    }
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::SetSpringConstant(double springConstant)
{
    mSpringConstant = springConstant;
}

template<unsigned DIM>
double ImmersedBoundaryMembraneElasticityForce<DIM>::GetSpringConstant()
{
    return mSpringConstant;
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::SetRestLengthMultiplier(double restLengthMultiplier)
{
    mRestLengthMultiplier = restLengthMultiplier;
}

template<unsigned DIM>
double ImmersedBoundaryMembraneElasticityForce<DIM>::GetRestLengthMultiplier()
{
    return mRestLengthMultiplier;
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SpringConstant>" << mSpringConstant << "</SpringConstant>\n";
    *rParamsFile << "\t\t\t<RestLengthMultiplier>" << mRestLengthMultiplier << "</RestLengthMultiplier>\n";
    *rParamsFile << "\t\t\t<BasementSpringConstantModifier>" << mBasementSpringConstantModifier << "</BasementSpringConstantModifier>\n";
    *rParamsFile << "\t\t\t<BasementRestLengthModifier>" << mBasementRestLengthModifier << "</BasementRestLengthModifier>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

// Explicit instantiation
template class ImmersedBoundaryMembraneElasticityForce<1>;
template class ImmersedBoundaryMembraneElasticityForce<2>;
template class ImmersedBoundaryMembraneElasticityForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryMembraneElasticityForce)
