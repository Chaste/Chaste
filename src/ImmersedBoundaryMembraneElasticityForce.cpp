/*

Copyright (c) 2005-2015, University of Oxford.
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
#include "Debug.hpp"

template<unsigned DIM>
ImmersedBoundaryMembraneElasticityForce<DIM>::ImmersedBoundaryMembraneElasticityForce(ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
        : AbstractImmersedBoundaryForce<DIM>(),
          mpCellPopulation(&rCellPopulation),
          mpMesh(&(rCellPopulation.rGetMesh())),
          mSpringConst(1e8),
          mRestLength(0.25 * mpMesh->GetCharacteristicNodeSpacing()),
          mBasementSpringConstantModifier(2.0),
          mBasementRestLengthModifier(0.5)
{
    // First verify that all elements have the same number of attributes
    mReferenceLocationInAttributesVector = mpMesh->GetElement(0)->GetNumElementAttributes();
    for (unsigned elem_idx = 1 ; elem_idx < mpMesh->GetNumElements() ; elem_idx++)
    {
        if (mReferenceLocationInAttributesVector != mpMesh->GetElement(elem_idx)->GetNumElementAttributes())
        {
            EXCEPTION("All elements must have the same number of attributes to use this force class.");
        }
    }

    /*
     * We split the nodes into three categories: basal, apical, and lateral.  We keep this information in the attribute
     * called region, with 0, 1, and 2 representing basal, apical, and lateral respectively.
     */
    TagNodeRegions();

    /*
     * We keep track of the initial size of the apical and basal sides.  This will be the initial distance between the
     * corners, which are storred by the element.
     *
     * Corners are represented as follows, and stored as four consecutive element attributes:
     *
     *     Apical
     *     0-----1
     *     |     |
     *     |     |
     *     |     |
     *     |     |
     *     |     |
     *     3-----2
     *      Basal
     *
     * The two element attributes store the starting distance between the apical corners and basal corners, giving
     * us:
     *
     * Attribute i:   Initial distance between apical corners
     *           i+1: Initial distance between basal corners
     */
    TagApicalAndBasalLengths();
}

template<unsigned DIM>
ImmersedBoundaryMembraneElasticityForce<DIM>::ImmersedBoundaryMembraneElasticityForce()
{
}

template<unsigned DIM>
ImmersedBoundaryMembraneElasticityForce<DIM>::~ImmersedBoundaryMembraneElasticityForce()
{
}

template<unsigned DIM>
double ImmersedBoundaryMembraneElasticityForce<DIM>::GetApicalLengthForElement(unsigned elemIndex)
{
    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mReferenceLocationInAttributesVector;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location];
}

template<unsigned DIM>
double ImmersedBoundaryMembraneElasticityForce<DIM>::GetBasalLengthForElement(unsigned elemIndex)
{
    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mReferenceLocationInAttributesVector + 1;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location];
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::AddForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = mpMesh->GetElementIteratorBegin();
         elem_it != mpMesh->GetElementIteratorEnd();
         ++elem_it)
    {
        // Get index and number of nodes of current element
        unsigned elem_idx = elem_it->GetIndex();
        unsigned num_nodes = elem_it->GetNumNodes();

        // Helper variables
        double normed_dist;
        c_vector<double, DIM> aggregate_force;

        // Make a vector to store the force on node i+1 from node i
        std::vector<c_vector<double, DIM> > elastic_force_to_next_node(num_nodes);

        double spring_constant = mSpringConst;
        double rest_length = mRestLength;

        /*
         * Here we make any necessary modifications to the spring properties
         */

        // The basement lamina, if present, will have different properties
        if (elem_idx == mpMesh->GetMembraneIndex())
        {
            spring_constant *= mBasementSpringConstantModifier;
            rest_length *= mBasementRestLengthModifier;
        }

        // Loop over nodes and calculate the force exerted on node i+1 by node i
        for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
        {
            // Index of the next node, calculated modulo number of nodes in this element
            unsigned next_idx = (node_idx + 1) % num_nodes;

            double modified_spring_constant = spring_constant;
            double modified_rest_length = rest_length;

            // Hooke's law linear spring force
            elastic_force_to_next_node[node_idx] = mpMesh->GetVectorFromAtoB(elem_it->GetNodeLocation(next_idx), elem_it->GetNodeLocation(node_idx));
            normed_dist = norm_2(elastic_force_to_next_node[node_idx]);
            elastic_force_to_next_node[node_idx] *= modified_spring_constant * (normed_dist - modified_rest_length) / normed_dist;
        }

        // Add the contributions of springs adjacent to each node
        for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
        {
            // Index of previous node, but -1%n doesn't work, so we add num_nodes when calculating
            unsigned prev_idx = (node_idx + num_nodes - 1) % num_nodes;

            aggregate_force = elastic_force_to_next_node[prev_idx] - elastic_force_to_next_node[node_idx];

            // Add the aggregate force contribution to the node
            elem_it->GetNode(node_idx)->AddAppliedForceContribution(aggregate_force);
        }

        // Add force contributions from apical and basal surfaces
        if (elem_idx != mpMesh->GetMembraneIndex())
        {
            std::vector<Node<DIM>*> r_corners = elem_it->rGetCornerNodes();

            // Apical surface
            c_vector<double, DIM> apical_force = mpMesh->GetVectorFromAtoB(r_corners[0]->rGetLocation(),
                                                                           r_corners[1]->rGetLocation());
            normed_dist = norm_2(apical_force);
            apical_force *=
                    1.0 * mSpringConst * (normed_dist - GetApicalLengthForElement(elem_idx)) / normed_dist;

            r_corners[0]->AddAppliedForceContribution(apical_force);
            apical_force *= -1.0;
            r_corners[1]->AddAppliedForceContribution(apical_force);

            // Basal surface
            c_vector<double, DIM> basal_force = mpMesh->GetVectorFromAtoB(r_corners[3]->rGetLocation(),
                                                                          r_corners[2]->rGetLocation());

            normed_dist = norm_2(basal_force);
            basal_force *=
                    1.0 * mSpringConst * (normed_dist - GetBasalLengthForElement(elem_idx)) / normed_dist;

            r_corners[3]->AddAppliedForceContribution(basal_force);
            basal_force *= -1.0;
            r_corners[2]->AddAppliedForceContribution(basal_force);
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::TagNodeRegions()
{
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = mpMesh->GetElementIteratorBegin();
         elem_it != mpMesh->GetElementIteratorEnd();
         ++elem_it)
    {
        // Basement lamina nodes are all basal
        if (mpMesh->GetMembraneIndex() == elem_it->GetIndex())
        {
            for (unsigned node_idx = 0 ; node_idx < elem_it->GetNumNodes() ; node_idx++)
            {
                elem_it->GetNode(node_idx)->SetRegion(msBas);
            }
        }
        else // not the basal lamina
        {
            // Nodes are ordered anti-clockwise, so nodes to the first corner will be lateral, then apical, lateral,
            // basal, lateral

            std::vector<Node<DIM>*> r_corners = elem_it->rGetCornerNodes();

            unsigned change_1 = elem_it->GetNodeLocalIndex(r_corners[1]->GetIndex());
            unsigned change_2 = elem_it->GetNodeLocalIndex(r_corners[0]->GetIndex()) + 1;
            unsigned change_3 = elem_it->GetNodeLocalIndex(r_corners[3]->GetIndex());
            unsigned change_4 = elem_it->GetNodeLocalIndex(r_corners[2]->GetIndex()) + 1;

            for (unsigned node_idx = 0 ; node_idx < change_1 ; node_idx++)
            {
                elem_it->GetNode(node_idx)->SetRegion(msLat);
            }
            for (unsigned node_idx = change_1 ; node_idx < change_2 ; node_idx++)
            {
                elem_it->GetNode(node_idx)->SetRegion(msApi);
            }
            for (unsigned node_idx = change_2 ; node_idx < change_3 ; node_idx++)
            {
                elem_it->GetNode(node_idx)->SetRegion(msLat);
            }
            for (unsigned node_idx = change_3 ; node_idx < change_4 ; node_idx++)
            {
                elem_it->GetNode(node_idx)->SetRegion(msBas);
            }
            for (unsigned node_idx = change_4 ; node_idx < elem_it->GetNumNodes() ; node_idx++)
            {
                elem_it->GetNode(node_idx)->SetRegion(msLat);
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::TagApicalAndBasalLengths()
{
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = mpMesh->GetElementIteratorBegin();
         elem_it != mpMesh->GetElementIteratorEnd();
         ++elem_it)
    {
        if (elem_it->GetIndex() == mpMesh->GetMembraneIndex())
        {
            // These lengths are irrelevant, but we will add them as zero to maintain the same number of element
            // attributes throughout
            elem_it->AddElementAttribute(0.0); //apical-length
            elem_it->AddElementAttribute(0.0); //basal-length
        }
        else
        {
            // Elements start roughly rectangular, so the correct apical and basal lengths are the width of the element

            unsigned half_way = unsigned(elem_it->GetNumNodes() / 2.0);

            double elem_width = fabs( mpMesh->GetVectorFromAtoB(elem_it->GetNode(0)->rGetLocation(),
                                                                elem_it->GetNode(half_way)->rGetLocation())[0] );


            elem_it->AddElementAttribute(elem_width);
            elem_it->AddElementAttribute(elem_width);
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
//    *rParamsFile << "\t\t\t<RestLength>" << mRestLength << "</RestLength>\n";
//    *rParamsFile << "\t\t\t<SpringConstant>" << mSpringConstant << "</SpringConstant>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ImmersedBoundaryMembraneElasticityForce<1>;
template class ImmersedBoundaryMembraneElasticityForce<2>;
template class ImmersedBoundaryMembraneElasticityForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryMembraneElasticityForce)
