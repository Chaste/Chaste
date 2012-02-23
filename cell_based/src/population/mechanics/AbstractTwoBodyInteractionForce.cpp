/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "AbstractTwoBodyInteractionForce.hpp"

template<unsigned DIM>
AbstractTwoBodyInteractionForce<DIM>::AbstractTwoBodyInteractionForce()
   : AbstractForce<DIM>(),
     mUseCutOffLength(false),
     mMechanicsCutOffLength(DBL_MAX)
{
}

template<unsigned DIM>
bool AbstractTwoBodyInteractionForce<DIM>::GetUseCutOffLength()
{
    return mUseCutOffLength;
}

template<unsigned DIM>
void AbstractTwoBodyInteractionForce<DIM>::SetCutOffLength(double cutOffLength)
{
    assert(cutOffLength > 0.0);
    mUseCutOffLength = true;
    mMechanicsCutOffLength = cutOffLength;
}

template<unsigned DIM>
double AbstractTwoBodyInteractionForce<DIM>::GetCutOffLength()
{
    return mMechanicsCutOffLength;
}

template<unsigned DIM>
void AbstractTwoBodyInteractionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                                AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a subclass of AbstractCentreBasedCellPopulation
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("Subclasses of AbstractTwoBodyInteractionForce are to be used with subclasses of AbstractCentreBasedCellPopulation only");
    }

    if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        MeshBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

        // Iterate over all springs and add force contributions
        for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
             spring_iterator != p_static_cast_cell_population->SpringsEnd();
             ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

            // Calculate the force between nodes
            c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation);

            // Add the force contribution to each node
            rForces[nodeB_global_index] -= force;
            rForces[nodeA_global_index] += force;
        }
    }
    else
    {
        std::set< std::pair<Node<DIM>*, Node<DIM>* > >& r_node_pairs = (static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))->rGetNodePairs();

//        assert(DIM==2); // 3d boxes not implemented yet - if fails nightly uncomment the double node loop below
                        // and use that for the 3d case
        for (typename std::set< std::pair<Node<DIM>*, Node<DIM>* > >::iterator iter = r_node_pairs.begin();
            iter != r_node_pairs.end();
            iter++)
        {
            std::pair<Node<DIM>*, Node<DIM>* > pair = *iter;

            unsigned node_a_index = pair.first->GetIndex();
            unsigned node_b_index = pair.second->GetIndex();

            // Calculate the force between nodes
            c_vector<double, DIM> force = CalculateForceBetweenNodes(node_a_index, node_b_index, rCellPopulation);
            for (unsigned j=0; j<DIM; j++)
            {
                assert(!std::isnan(force[j]));
            }

            // Add the force contribution to each node
            rForces[node_a_index] += force;
            rForces[node_b_index] -= force;
        }

//        // Iterate over nodes
//        for (unsigned node_a_index=0; node_a_index<rCellPopulation.GetNumNodes(); node_a_index++)
//        {
//            // Iterate over nodes
//            for (unsigned node_b_index=node_a_index+1; node_b_index<rCellPopulation.GetNumNodes(); node_b_index++)
//            {
//                // Calculate the force between nodes
//                c_vector<double, DIM> force = CalculateForceBetweenNodes(node_a_index, node_b_index, rCellPopulation);
//                for (unsigned j=0; j<DIM; j++)
//                {
//                    assert(!std::isnan(force[j]));
//                }
//
//                // Add the force contribution to each node
//                rForces[node_a_index] += force;
//                rForces[node_b_index] -= force;
//            }
//        }
    }
}

template<unsigned DIM>
void AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<UseCutOffLength>" << mUseCutOffLength << "</UseCutOffLength>\n";
    *rParamsFile << "\t\t\t<CutOffLength>" << mMechanicsCutOffLength << "</CutOffLength>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AbstractTwoBodyInteractionForce<1>;
template class AbstractTwoBodyInteractionForce<2>;
template class AbstractTwoBodyInteractionForce<3>;
