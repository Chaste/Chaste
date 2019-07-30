/*

Copyright (c) 2005-2019, University of Oxford.
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::AbstractTwoBodyInteractionForce()
   : AbstractForce<ELEMENT_DIM,SPACE_DIM>(),
     mUseCutOffLength(false),
     mMechanicsCutOffLength(DBL_MAX)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::GetUseCutOffLength()
{
    return mUseCutOffLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::SetCutOffLength(double cutOffLength)
{
    assert(cutOffLength > 0.0);
    mUseCutOffLength = true;
    mMechanicsCutOffLength = cutOffLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::GetCutOffLength()
{
    return mMechanicsCutOffLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // Throw an exception message if not using a subclass of AbstractCentreBasedCellPopulation
    if (dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("Subclasses of AbstractTwoBodyInteractionForce are to be used with subclasses of AbstractCentreBasedCellPopulation only");
    }

    ///\todo this could be tidied by using the rGetNodePairs for all populations and moving the below calculation into the MutableMesh.
    if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
    {
        MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

        // Iterate over all springs and add force contributions
        for (typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
             spring_iterator != p_static_cast_cell_population->SpringsEnd();
             ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

            // Calculate the force between nodes
            c_vector<double, SPACE_DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation);

            // Add the force contribution to each node
            c_vector<double, SPACE_DIM> negative_force = -1.0*force;
            spring_iterator.GetNodeB()->AddAppliedForceContribution(negative_force);
            spring_iterator.GetNodeA()->AddAppliedForceContribution(force);
        }
    }
    else    // This is a NodeBasedCellPopulation
    {
        AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

        std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >& r_node_pairs = p_static_cast_cell_population->rGetNodePairs();

        for (typename std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >::iterator iter = r_node_pairs.begin();
            iter != r_node_pairs.end();
            iter++)
        {
            std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > pair = *iter;

            unsigned node_a_index = pair.first->GetIndex();
            unsigned node_b_index = pair.second->GetIndex();

            // Calculate the force between nodes
            c_vector<double, SPACE_DIM> force = CalculateForceBetweenNodes(node_a_index, node_b_index, rCellPopulation);
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                assert(!std::isnan(force[j]));
            }

            // Add the force contribution to each node
            c_vector<double, SPACE_DIM> negative_force = -1.0*force;
            pair.first->AddAppliedForceContribution(force);
            pair.second->AddAppliedForceContribution(negative_force);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<UseCutOffLength>" << mUseCutOffLength << "</UseCutOffLength>\n";
    *rParamsFile << "\t\t\t<CutOffLength>" << mMechanicsCutOffLength << "</CutOffLength>\n";

    // Call method on direct parent class
    AbstractForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile)
{
    *pVizSetupFile << "Cutoff\t" << mMechanicsCutOffLength << "\n";
}

// Explicit instantiation
template class AbstractTwoBodyInteractionForce<1,1>;
template class AbstractTwoBodyInteractionForce<1,2>;
template class AbstractTwoBodyInteractionForce<2,2>;
template class AbstractTwoBodyInteractionForce<1,3>;
template class AbstractTwoBodyInteractionForce<2,3>;
template class AbstractTwoBodyInteractionForce<3,3>;
