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

#include "RepulsionForce.hpp"

template<unsigned DIM>
RepulsionForce<DIM>::RepulsionForce()
   : GeneralisedLinearSpringForce<DIM>()
{
}

template<unsigned DIM>
void RepulsionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                               AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a NodeBasedCellPopulation
    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("RepulsionForce is to be used with a NodeBasedCellPopulation only");
    }

    std::set< std::pair<Node<DIM>*, Node<DIM>* > >& r_node_pairs = (static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))->rGetNodePairs();

    for (typename std::set< std::pair<Node<DIM>*, Node<DIM>* > >::iterator iter = r_node_pairs.begin();
        iter != r_node_pairs.end();
        iter++)
    {
        std::pair<Node<DIM>*, Node<DIM>* > pair = *iter;

        unsigned node_a_index = pair.first->GetIndex();
        unsigned node_b_index = pair.second->GetIndex();

        // Get the node locations
        c_vector<double, DIM> node_a_location = rCellPopulation.GetNode(node_a_index)->rGetLocation();
        c_vector<double, DIM> node_b_location = rCellPopulation.GetNode(node_b_index)->rGetLocation();

        // Get the unit vector parallel to the line joining the two nodes
        c_vector<double, DIM> unit_difference;

        unit_difference = node_b_location - node_a_location;


        double rest_length = 1.0;
        if (norm_2(unit_difference) < rest_length)
        {
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
    }
}

template<unsigned DIM>
void RepulsionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call direct parent class
    GeneralisedLinearSpringForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class RepulsionForce<1>;
template class RepulsionForce<2>;
template class RepulsionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RepulsionForce)
