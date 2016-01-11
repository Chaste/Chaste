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

#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryElement.hpp"

template<unsigned DIM>
ImmersedBoundaryMembraneElasticityForce<DIM>::ImmersedBoundaryMembraneElasticityForce(ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
        : AbstractImmersedBoundaryForce<DIM>(),
          mpCellPopulation(&rCellPopulation),
          mpMesh(&(rCellPopulation.rGetMesh())),
mSpringConstant()
{
    // First verify that all elements have the same number of attributes
    unsigned num_elem_attributes = mpMesh->GetElement(0)->GetNumElementAttributes();
    for (unsigned elem_idx = 0 ; elem_idx < mpMesh->GetNumElements() ; elem_idx++)
    {
        if (num_elem_attributes != mpMesh->GetElement(elem_idx)->GetNumElementAttributes())
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
     * We calculate the 'corners' of each element, in order to alter behaviour on apical, lateral, and basal regions
     * separately.
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
     */
//    TagElementCorners();
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
void ImmersedBoundaryMembraneElasticityForce<DIM>::AddForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    ImmersedBoundaryMesh<DIM,DIM>* p_mesh = &(mpCellPopulation->rGetMesh());

    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_iter = p_mesh->GetElementIteratorBegin();
         elem_iter != p_mesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get number of nodes in current element
        unsigned num_nodes = elem_iter->GetNumNodes();
        assert(num_nodes > 0);

        // Get spring parameters (owned by the elements as they may differ within the element population)
        double spring_constant = elem_iter->GetMembraneSpringConstant();
        double rest_length = elem_iter->GetMembraneRestLength();

        // Helper variables
        double normed_dist;
        c_vector<double, DIM> aggregate_force;

        // Make a vector to store the force on node i+1 from node i
        std::vector<c_vector<double, DIM> > elastic_force_to_next_node(num_nodes);

        // Loop over nodes and calculate the force exerted on node i+1 by node i
        for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
        {
            // Index of the next node, calculated modulo number of nodes in this element
            unsigned next_idx = (node_idx + 1) % num_nodes;

            double modified_spring_constant = spring_constant;
            double modified_rest_length = rest_length;

            // If the node is apical or basal, increase the spring constant
            if (elem_iter->GetNode(node_idx)->GetRegion() < 2)
            {
                modified_spring_constant *= 10.0;
                modified_rest_length *= 4.0;
            }

            // Hooke's law linear spring force
            elastic_force_to_next_node[node_idx] = p_mesh->GetVectorFromAtoB(elem_iter->GetNodeLocation(next_idx), elem_iter->GetNodeLocation(node_idx));
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
            elem_iter->GetNode(node_idx)->AddAppliedForceContribution(aggregate_force);
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::TagNodeRegions()
{
    ImmersedBoundaryMesh<DIM,DIM>* p_mesh = &(mpCellPopulation->rGetMesh());

    for (unsigned elem_idx = 0 ; elem_idx < p_mesh->GetNumElements() ; elem_idx++)
    {
        ImmersedBoundaryElement<DIM,DIM>* p_this_elem = p_mesh->GetElement(elem_idx);

        // Basement lamina nodes are all basal
        if (p_mesh->GetMembraneIndex() == p_this_elem->GetIndex())
        {
            for (unsigned node_idx = 0 ; node_idx < p_this_elem->GetNumNodes() ; node_idx++)
            {
                p_this_elem->GetNode(node_idx)->SetRegion(2);
            }
        }
        else // not the basal lamina
        {
            /*
             * Because cells are initialised as roughly rectangular structures with equally spaced nodes, the correct
             * number of basal (or apical) nodes will be (roughly) 0.5 * num_nodes / (1 + aspect ratio).
             *
             * We identify which cells will be apical and basal by sorting the locations of each node and calculating
             * the correct threshold values, which we then compare against when assigning the region.
             */
            unsigned num_nodes = p_this_elem->GetNumNodes();
            double aspect_ratio = p_mesh->GetElongationShapeFactorOfElement(elem_idx);

            unsigned num_basal_nodes = std::floor(0.5 * ( double(num_nodes) / (1.0 + aspect_ratio) ));

            // Check we have more than one, and fewer than half, of the nodes to be marked as basal
            assert(num_basal_nodes > 1);
            assert(num_basal_nodes < std::floor(double(num_nodes) / 2.0) );

            std::vector<double> node_y_locations;

            for (unsigned node_idx = 0 ; node_idx < p_this_elem->GetNumNodes() ; node_idx++)
            {
                node_y_locations.push_back(p_this_elem->GetNode(node_idx)->rGetLocation()[1]);
            }

            std::sort(node_y_locations.begin(), node_y_locations.end());

            double low_threshold  = 0.5 * (node_y_locations[num_basal_nodes - 1] + node_y_locations[num_basal_nodes]);
            double high_threshold = 0.5 * (node_y_locations[num_nodes - num_basal_nodes] + node_y_locations[num_nodes - num_basal_nodes - 1]);

            assert(low_threshold < high_threshold);

            for (unsigned node_idx = 0 ; node_idx < p_this_elem->GetNumNodes() ; node_idx++)
            {
                double node_y_location = p_this_elem->GetNode(node_idx)->rGetLocation()[1];

                if (node_y_location < low_threshold)
                {
                    // Node will be basal (region 0)
                    p_this_elem->GetNode(node_idx)->SetRegion(0);
                }
                else if (node_y_location > high_threshold)
                {
                    // Node will be apical (region 1)
                    p_this_elem->GetNode(node_idx)->SetRegion(1);
                }
                else
                {
                    // Node will be lateral (region 2)
                    p_this_elem->GetNode(node_idx)->SetRegion(2);
                }
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryMembraneElasticityForce<DIM>::TagElementCorners()
{
//    ImmersedBoundaryMesh<DIM,DIM>* p_mesh = &(mpCellPopulation->rGetMesh());
//
//    // First loop through all elements to check they have the same number of attributes
//    assert(p_mesh->GetNumElements() > 0);
//    unsigned num_elem_attributes = p_mesh->GetElement(0)->GetNumElementAttributes();
//
//    for (unsigned elem_idx = 1 ; elem_idx < p_mesh->GetNumElements() ; elem_idx++)
//    {
//        if (p_mesh->GetElement(elem_idx)->GetNumElementAttributes() != num_elem_attributes)
//        {
//            EXCEPTION("This class requires each element to have the same number of attributes");
//        }
//    }
//
//    // Set up corner locations in the attribute vector, so AddForceContribution knows where to look
//    mCornerLocationsInAttributeVector.push_back(num_elem_attributes);
//    mCornerLocationsInAttributeVector.push_back(num_elem_attributes + 1);
//    mCornerLocationsInAttributeVector.push_back(num_elem_attributes + 2);
//    mCornerLocationsInAttributeVector.push_back(num_elem_attributes + 3);
//
//    // Set up rest-length locations in the attribute vector, so AddForceContribution knows where to look
//    mRestLengthLocationsInAttributeVector.push_back(num_elem_attributes + 4);
//    mRestLengthLocationsInAttributeVector.push_back(num_elem_attributes + 5);
//
//    /*
//     * Loop through elements and set corner locations.
//     *
//     * Corner 0 will be the left-most  apical node.
//     * Corner 1 will be the right-most apical node.
//     * Corner 2 will be the right-most basal node.
//     * Corner 3 will be the left-most  basal node.
//     */
//    for (unsigned elem_idx = 0 ; elem_idx < p_mesh->GetNumElements() ; elem_idx++)
//    {
//        ImmersedBoundaryElement<DIM,DIM>* p_this_elem = p_mesh->GetElement(elem_idx);
//
//        if (p_mesh->GetMembraneIndex() == p_this_elem->GetIndex())
//        {
//            /*
//             * The basement lamina need not have any attributes associated with it, as it is treated differently when
//             * calculating forces.  However, we require each element to have the same number of attributes, so we add
//             * arbitrary attributes (with no meaning) to the vector.
//             */
//
//            // Basement lamina will have no corners, so we set each corner to 0
//            p_this_elem->AddElementAttribute(0.0);
//            p_this_elem->AddElementAttribute(0.0);
//            p_this_elem->AddElementAttribute(0.0);
//            p_this_elem->AddElementAttribute(0.0);
//
//            // Basement lamina will have rest-length attributes set to 1.0
//            p_this_elem->AddElementAttribute(1.0);
//            p_this_elem->AddElementAttribute(1.0);
//        }
//        else // not the basal lamina
//        {
//            double left_most_apical  = DBL_MAX;
//            double right_most_apical = -DBL_MAX;
//            double left_most_basal   = DBL_MAX;
//            double right_most_basal  = -DBL_MAX;
//
//            unsigned left_apical_idx;
//            unsigned right_apical_idx;
//            unsigned left_basal_idx;
//            unsigned right_basal_idx;
//
//            c_vector<double, DIM> origin = zero_vector<double>(DIM);
//
//            // Calculate the node indices corresponding to the left-most and right-most apical and basal nodes in the
//            // current element
//            for (unsigned node_idx = 0 ; node_idx < p_this_elem->GetNumNodes() ; node_idx++)
//            {
//                Node<DIM>* p_this_node = p_this_elem->GetNode(node_idx);
//
//                if (p_this_node->GetRegion() == 1) //apical
//                {
//                    // Get the deviation of the current location from the origin, taking into account periodicity
//                    double x_deviation = p_mesh->GetVectorFromAtoB(origin, p_this_elem->GetNode(node_idx)->rGetLocation())[0];
//
//                    if(x_deviation < left_most_apical)
//                    {
//                        left_most_apical = x_deviation;
//                        left_apical_idx = node_idx;
//                    }
//                    if(x_deviation > right_most_apical)
//                    {
//                        right_most_apical = x_deviation;
//                        right_apical_idx = node_idx;
//                    }
//                }
//                else if (p_this_node->GetRegion() == 0) //basal
//                {
//                    // Get the deviation of the current location from the origin, taking into account periodicity
//                    double x_deviation = p_mesh->GetVectorFromAtoB(origin, p_this_elem->GetNode(node_idx)->rGetLocation())[0];
//
//                    if(x_deviation < left_most_basal)
//                    {
//                        left_most_basal = x_deviation;
//                        left_basal_idx = node_idx;
//                    }
//                    if(x_deviation > right_most_basal)
//                    {
//                        right_most_basal = x_deviation;
//                        right_basal_idx = node_idx;
//                    }
//                }
//                else
//                {
//                    // Lateral - nothing to do
//                }
//            }
//
//            double apical_rest_length = norm_2(p_mesh->GetVectorFromAtoB(p_this_elem->GetNode(left_apical_idx)->rGetLocation()),
//                                               p_mesh->GetVectorFromAtoB(p_this_elem->GetNode(right_apical_idx)->rGetLocation()));
//
//            // Basement lamina will have no corners, so we set each corner to 0
//            p_this_elem->AddElementAttribute(left_apical_idx);
//            p_this_elem->AddElementAttribute(0.0);
//            p_this_elem->AddElementAttribute(0.0);
//            p_this_elem->AddElementAttribute(0.0);
//
//            // Basement lamina will have rest-length attributes set to 1.0
//            p_this_elem->AddElementAttribute(1.0);
//            p_this_elem->AddElementAttribute(1.0);
//        }
//    }
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
