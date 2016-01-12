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
          mSpringConst(1e9),
          mRestLength(0.25 * mpMesh->GetCharacteristicNodeSpacing()),
          mBasementSpringConstantModifier(2.0),
          mBasementRestLengthModifier(0.5)
{
    // First verify that all elements have the same number of attributes
    mCurrentLocationInElementAttributesVector = mpMesh->GetElement(0)->GetNumElementAttributes();
    for (unsigned elem_idx = 1 ; elem_idx < mpMesh->GetNumElements() ; elem_idx++)
    {
        if (mCurrentLocationInElementAttributesVector != mpMesh->GetElement(elem_idx)->GetNumElementAttributes())
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
     * The next two element attributes store the starting distance between the apical corners and basal corners, giving
     * us:
     *
     * Attribute i:   Left-apical-corner node index
     *           i+1: Right-apical-corner node index
     *           i+2: Right-basal-corner node index
     *           i+3: Left-basal-corner node index
     *           i+4: Initial distance between apical corners
     *           i+5: Initial distance between basal corners
     */
    TagElementCorners();
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
unsigned ImmersedBoundaryMembraneElasticityForce<DIM>::GetLeftApicalCornerNodeIndexForElement(unsigned elemIndex)
{
    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mCurrentLocationInElementAttributesVector;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return unsigned(mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location]);
}

template<unsigned DIM>
unsigned ImmersedBoundaryMembraneElasticityForce<DIM>::GetRightApicalCornerNodeIndexForElement(unsigned elemIndex)
{
    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mCurrentLocationInElementAttributesVector + 1;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return unsigned(mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location]);
}

template<unsigned DIM>
unsigned ImmersedBoundaryMembraneElasticityForce<DIM>::GetRightBasalCornerNodeIndexForElement(unsigned elemIndex)
{
    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mCurrentLocationInElementAttributesVector + 2;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return unsigned(mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location]);
}

template<unsigned DIM>
unsigned ImmersedBoundaryMembraneElasticityForce<DIM>::GetLeftBasalCornerNodeIndexForElement(unsigned elemIndex)
{
    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mCurrentLocationInElementAttributesVector + 3;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return unsigned(mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location]);
}

template<unsigned DIM>
double ImmersedBoundaryMembraneElasticityForce<DIM>::GetApicalLengthForElement(unsigned elemIndex)
{
    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mCurrentLocationInElementAttributesVector + 4;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location];
}

template<unsigned DIM>
double ImmersedBoundaryMembraneElasticityForce<DIM>::GetBasalLengthForElement(unsigned elemIndex)
{
    // Calculate correct location in attributes vector and check it's valid
    unsigned attribute_location = mCurrentLocationInElementAttributesVector + 5;
    assert(attribute_location < mpMesh->GetElement(elemIndex)->GetNumElementAttributes());

    return mpMesh->GetElement(elemIndex)->rGetElementAttributes()[attribute_location];
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
        if (elem_iter->GetIndex() == mpMesh->GetMembraneIndex())
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

        // Add force contributions from apical and basal parts
        if (elem_iter->GetIndex() != mpMesh->GetMembraneIndex())
        {
            // Apical nodes
            Node<DIM> *p_apical_left = elem_iter->GetNode(GetLeftApicalCornerNodeIndexForElement(elem_iter->GetIndex()));
            Node<DIM> *p_apical_right = elem_iter->GetNode(GetRightApicalCornerNodeIndexForElement(elem_iter->GetIndex()));

            c_vector<double, DIM> apical_force = p_mesh->GetVectorFromAtoB(p_apical_left->rGetLocation(),
                                                                           p_apical_right->rGetLocation());
            normed_dist = norm_2(apical_force);
            apical_force *=
                    0.1 * mSpringConst * (normed_dist - GetApicalLengthForElement(elem_iter->GetIndex())) / normed_dist;

            p_apical_left->AddAppliedForceContribution(apical_force);
            apical_force *= -1.0;
            p_apical_right->AddAppliedForceContribution(apical_force);


            // Basal nodes
            Node<DIM> *p_basal_left = elem_iter->GetNode(GetLeftBasalCornerNodeIndexForElement(elem_iter->GetIndex()));
            Node<DIM> *p_basal_right = elem_iter->GetNode(GetRightBasalCornerNodeIndexForElement(elem_iter->GetIndex()));

            c_vector<double, DIM> basal_force = p_mesh->GetVectorFromAtoB(p_basal_left->rGetLocation(),
                                                                          p_basal_right->rGetLocation());
            PRINT_VECTOR(basal_force);
            PRINT_VECTOR(p_basal_left->rGetLocation());
            PRINT_VECTOR(p_basal_right->rGetLocation());
            normed_dist = norm_2(basal_force);
            basal_force *=
                    0.0 * mSpringConst * (normed_dist - GetBasalLengthForElement(elem_iter->GetIndex())) / normed_dist;

            p_basal_left->AddAppliedForceContribution(basal_force);
            basal_force *= -1.0;
            p_basal_right->AddAppliedForceContribution(basal_force);
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
    /*
     * Loop through elements and set corner locations.
     *
     * Corner 0 will be the left-most  apical node.
     * Corner 1 will be the right-most apical node.
     * Corner 2 will be the right-most basal node.
     * Corner 3 will be the left-most  basal node.
     */
    for (unsigned elem_idx = 0 ; elem_idx < mpMesh->GetNumElements() ; elem_idx++)
    {
        if (elem_idx == mpMesh->GetMembraneIndex())
        {
            // Corners are irrelevant for the basement lamina, so set all six attributes to zero
            mpMesh->GetElement(elem_idx)->AddElementAttribute(0.0); //L-A
            mpMesh->GetElement(elem_idx)->AddElementAttribute(0.0); //R-A
            mpMesh->GetElement(elem_idx)->AddElementAttribute(0.0); //R-B
            mpMesh->GetElement(elem_idx)->AddElementAttribute(0.0); //L-B
            mpMesh->GetElement(elem_idx)->AddElementAttribute(0.0); //A-length
            mpMesh->GetElement(elem_idx)->AddElementAttribute(0.0); //B-length
        }
        else
        {
            /*
             * Due to the pay elements are created using the palisade mesh generator, node 0 should be lateral, and
             * proceeding anticlockwise we should first come to a contiguous set of apical nodes, then lateral, then
             * basal.
             */

            if (mpMesh->GetElement(elem_idx)->GetNode(0)->GetRegion() != 2)
            {
                EXCEPTION("This class is intended only for use with the ImmersedBoundaryPalisadeMeshGenerator class");
            }

            unsigned apical_left = UINT_MAX;
            unsigned apical_right = UINT_MAX;
            unsigned basal_right = UINT_MAX;
            unsigned basal_left = UINT_MAX;

            unsigned num_nodes = mpMesh->GetElement(elem_idx)->GetNumNodes();

            unsigned last_region = 2; //lateral

            for (unsigned node_idx = 1 ; node_idx < num_nodes ; node_idx++)
            {
                unsigned this_region = mpMesh->GetElement(elem_idx)->GetNode(node_idx)->GetRegion();

                if (this_region == 1 && last_region == 2)  //lateral -> apical transition
                {
                    apical_right = node_idx;
                }
                else if (this_region == 2 && last_region == 1) // apical -> lateral transition
                {
                    apical_left = node_idx - 1;
                }
                else if (this_region == 0 && last_region == 2) // lateral -> basal transition
                {
                    basal_left = node_idx;
                }
                else if (this_region == 2 && last_region == 0) // basal -> lateral transition
                {
                    basal_right = node_idx - 1;
                }

                last_region = this_region;
            }

            if (apical_left == UINT_MAX || apical_right == UINT_MAX || basal_right == UINT_MAX || basal_left == UINT_MAX)
            {
                EXCEPTION("Some corner nodes not tagged.");
            }

            if (apical_right >= apical_left || apical_left >= basal_left || basal_left >= basal_right)
            {
                EXCEPTION("Corner nodes have not been tagged correctly.");
            }

            mpMesh->GetElement(elem_idx)->AddElementAttribute(double(apical_left));
            mpMesh->GetElement(elem_idx)->AddElementAttribute(double(apical_right));
            mpMesh->GetElement(elem_idx)->AddElementAttribute(double(basal_right));
            mpMesh->GetElement(elem_idx)->AddElementAttribute(double(basal_left));

            double apical_length = norm_2(mpMesh->GetElement(elem_idx)->GetNode(apical_left)->rGetLocation() -
                                          mpMesh->GetElement(elem_idx)->GetNode(apical_right)->rGetLocation());

            double basal_length  = norm_2(mpMesh->GetElement(elem_idx)->GetNode(basal_left)->rGetLocation() -
                                          mpMesh->GetElement(elem_idx)->GetNode(basal_right)->rGetLocation());

            mpMesh->GetElement(elem_idx)->AddElementAttribute(apical_length);
            mpMesh->GetElement(elem_idx)->AddElementAttribute(basal_length);
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
