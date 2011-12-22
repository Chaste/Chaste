/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
