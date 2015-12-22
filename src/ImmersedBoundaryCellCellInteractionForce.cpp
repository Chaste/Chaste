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

#include "ImmersedBoundaryCellCellInteractionForce.hpp"

#include "ImmersedBoundaryElement.hpp"
#include "Debug.hpp"

template<unsigned DIM>
ImmersedBoundaryCellCellInteractionForce<DIM>::ImmersedBoundaryCellCellInteractionForce(ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
        : AbstractImmersedBoundaryForce<DIM>(),
          mpCellPopulation(&rCellPopulation),
          mpMesh(&(rCellPopulation.rGetMesh())),
          mSpringConst(1e3),
          mRestLength(0.25 * rCellPopulation.GetInteractionDistance())
{
    /*
     * This force class calculates the force between pairs of nodes in different immersed boundaries.  Each node must
     * therefore store a dimensionless parameter representing the quantity of different transmembrane proteins at that
     * location.  We attach these quantities as node attributes, and keep track of where in the node attributes vector
     * each protein concentration is stored.
     */

    // First verify that all nodes have the same number of attributes
    unsigned num_node_attributes = mpMesh->GetNode(0)->GetNumNodeAttributes();
    for (unsigned node_idx = 0 ; node_idx < mpMesh->GetNumNodes() ; node_idx++ )
    {
        if (num_node_attributes != mpMesh->GetNode(node_idx)->GetNumNodeAttributes())
        {
            EXCEPTION("All nodes must have teh same number of attributes to use this force class.");
        }
    }

    // Set up the number of proteins and keep track of where they will be stored in the node attributes vector
    mNumProteins = 3;
    for (unsigned protein_idx = 0 ; protein_idx < mNumProteins ; protein_idx++)
    {
        mProteinNodeAttributeLocations.push_back(num_node_attributes + protein_idx);
    }

    // Add protein attributes to each node
    for (unsigned node_idx = 0 ; node_idx < mpMesh->GetNumNodes() ; node_idx++)
    {
        for (unsigned protein_idx = 0; protein_idx < mNumProteins; protein_idx++)
        {
            mpMesh->GetNode(node_idx)->AddNodeAttribute(0.0);
        }
    }

    //Initialize protein levels
    InitializeProteinLevels();
}

template<unsigned DIM>
ImmersedBoundaryCellCellInteractionForce<DIM>::ImmersedBoundaryCellCellInteractionForce()
{
}

template<unsigned DIM>
ImmersedBoundaryCellCellInteractionForce<DIM>::~ImmersedBoundaryCellCellInteractionForce()
{
}

template<unsigned DIM>
void ImmersedBoundaryCellCellInteractionForce<DIM>::AddForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    UpdateProteinLevels();

    std::cout << "Timesteps elapsed: " << SimulationTime::Instance()->GetTimeStepsElapsed() << std::endl;

    // Helper variables for loop
    unsigned e_cad_idx = mProteinNodeAttributeLocations[0];
    unsigned p_cad_idx = mProteinNodeAttributeLocations[1];
    unsigned integrin_idx = mProteinNodeAttributeLocations[2];

    double normed_dist;
    double protein_mult;

    c_vector<double, DIM> vector_between_nodes;

    Node<DIM>* p_node_a;
    Node<DIM>* p_node_b;

    // Loop over all pairs of nodes that might be interacting
    for (unsigned pair = 0 ; pair < rNodePairs.size() ; pair++)
    {
        /*
         * Interactions only occur between different cells.  Since each node is only ever in a single cell, we can test
         * equality of the first ContainingElement.
         */
        if ( *(rNodePairs[pair].first->ContainingElementsBegin()) !=
             *(rNodePairs[pair].second->ContainingElementsBegin()) )
        {
            p_node_a = rNodePairs[pair].first;
            p_node_b = rNodePairs[pair].second;

            std::vector<double>& r_a_attribs = p_node_a->rGetNodeAttributes();
            std::vector<double>& r_b_attribs = p_node_b->rGetNodeAttributes();

            vector_between_nodes = mpMesh->GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());
            normed_dist = norm_2(vector_between_nodes);

            if (normed_dist < mpCellPopulation->GetInteractionDistance())
            {
                // The protein multiplier is a function of the levels of each protein in the current and comparison nodes
                protein_mult = std::min(r_a_attribs[e_cad_idx], r_b_attribs[e_cad_idx]) +
                               std::min(r_a_attribs[p_cad_idx], r_b_attribs[p_cad_idx]) +
                               std::max(r_a_attribs[integrin_idx], r_b_attribs[integrin_idx]);

                vector_between_nodes *= mSpringConst * protein_mult * (normed_dist - mRestLength) / normed_dist;
                p_node_a->AddAppliedForceContribution(vector_between_nodes);

                vector_between_nodes *= -1.0;
                p_node_b->AddAppliedForceContribution(vector_between_nodes);
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryCellCellInteractionForce<DIM>::InitializeProteinLevels()
{
    /*
     * We are thinking of the following proteins:
     *  * 0: E-cadherin
     *  * 1: P-cadherin
     *  * 2: Integrins
     */
    for (unsigned elem_idx = 0 ; elem_idx < mpMesh->GetNumElements() ; elem_idx++)
    {
        double e_cad = 0.0;
        double p_cad = 0.0;
        double integrin = 0.0;

        if (mpMesh->GetElement(elem_idx) == mpMesh->GetMembraneElement())
        {
            integrin = 100.0;
        }
        else
        {
            e_cad = 1.0;
        }

        for (unsigned node_idx = 0 ; node_idx < mpMesh->GetElement(elem_idx)->GetNumNodes() ; node_idx++)
        {
            std::vector<double>& r_node_attributes = mpMesh->GetElement(elem_idx)->GetNode(node_idx)->rGetNodeAttributes();

            r_node_attributes[mProteinNodeAttributeLocations[0]] += e_cad;
            r_node_attributes[mProteinNodeAttributeLocations[1]] += p_cad;
            r_node_attributes[mProteinNodeAttributeLocations[2]] += integrin;
        }
    }
}

template<unsigned DIM>
void ImmersedBoundaryCellCellInteractionForce<DIM>::UpdateProteinLevels()
{

}

template<unsigned DIM>
void ImmersedBoundaryCellCellInteractionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
//    *rParamsFile << "\t\t\t<RestLength>" << mRestLength << "</RestLength>\n";
//    *rParamsFile << "\t\t\t<SpringConstant>" << mSpringConstant << "</SpringConstant>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ImmersedBoundaryCellCellInteractionForce<1>;
template class ImmersedBoundaryCellCellInteractionForce<2>;
template class ImmersedBoundaryCellCellInteractionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryCellCellInteractionForce)
