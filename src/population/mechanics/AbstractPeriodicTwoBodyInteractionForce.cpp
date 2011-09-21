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

#include "AbstractPeriodicTwoBodyInteractionForce.hpp"
#include "Exception.hpp"

template<unsigned DIM>
AbstractPeriodicTwoBodyInteractionForce<DIM>::AbstractPeriodicTwoBodyInteractionForce()
   : AbstractTwoBodyInteractionForce<DIM>(),
     mPeriodicDomainWidth(DOUBLE_UNSET),
     mPeriodicDomainDepth(DOUBLE_UNSET),
     mpExtendedMesh(NULL)
{
}

template<unsigned DIM>
AbstractPeriodicTwoBodyInteractionForce<DIM>::~AbstractPeriodicTwoBodyInteractionForce()
{
    delete mpExtendedMesh;
}

template<unsigned DIM>
void AbstractPeriodicTwoBodyInteractionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                                        AbstractCellPopulation<DIM>& rCellPopulation)
{
    // If the width of the periodic domain has not been specified, use the initial width of the cell population
    if (mPeriodicDomainWidth == DOUBLE_UNSET)
    {
        mPeriodicDomainWidth = rCellPopulation.GetWidth(0);
    }

    mExtendedMeshNodeIndexMap.clear();

    // Create a vector of nodes for use in constructing mpExtendedMesh
    unsigned num_cells = rCellPopulation.GetNumRealCells();
    std::vector<Node<DIM>*> extended_nodes(DIM*num_cells);

    // We iterate over all cells in the population
    unsigned count = 0;

    switch (DIM)
    {
        case 1:
        {
            break;
        }
        case 2:
        {
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
            {
                // First, create and store a copy of this real node and cell
                unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                // Create a copy of the node corresponding to this cell and store it
                Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
                extended_nodes[count] = p_real_node;

                // Compute the location of the image node corresponding to this node
                c_vector<double,DIM> image_node_location = real_node_location;
                if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
                {
                    image_node_location[0] -= mPeriodicDomainWidth;
                }
                else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
                {
                    image_node_location[0] += mPeriodicDomainWidth;
                }

                // Create a copy of the node corresponding to this cell, suitable translated, and store it
                Node<DIM>* p_image_node = new Node<DIM>(num_cells+count, image_node_location);
                extended_nodes[num_cells+count] = p_image_node;

                // Populate mExtendedMeshNodeIndexMap
                mExtendedMeshNodeIndexMap[count] = real_node_index;
                mExtendedMeshNodeIndexMap[num_cells+count] = real_node_index;

                count++;
            }

            // We now construct mpExtendedMesh using extended_nodes
            mpExtendedMesh = new MutableMesh<DIM,DIM>(extended_nodes);

        	// Now loop over the extended mesh and calculate the force acting on real nodes
        	// (using the edge iterator ensures that each edge is visited exactly once)
            for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator = mpExtendedMesh->EdgesBegin();
                 edge_iterator != mpExtendedMesh->EdgesEnd();
                 ++edge_iterator)
            {
                unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
                unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

                c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation);

                // Apply this force to any real nodes (i.e. nodes whose indices are less than num_real_nodes)
                if (nodeA_global_index < num_cells)
                {
                    unsigned real_node_index_A = mExtendedMeshNodeIndexMap[nodeA_global_index];
                    rForces[real_node_index_A] += force;
                }
                if (nodeB_global_index < num_cells)
                {
                    unsigned real_node_index_B = mExtendedMeshNodeIndexMap[nodeB_global_index];
                    rForces[real_node_index_B] -= force;
                }
            }

			break;
        }
        case 3:
        {
            // If the width of the periodic domain has not been specified, use the initial width of the cell population
            if (mPeriodicDomainDepth == DOUBLE_UNSET)
            {
            	mPeriodicDomainDepth = rCellPopulation.GetWidth(1) + 1.0;
            }

            // First, extend the mesh in the x-direction

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
            {
                // First, create and store a copy of this real node and cell
                unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                // Create a copy of the node corresponding to this cell and store it
                Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
                extended_nodes[count] = p_real_node;

                // Compute the location of the image node corresponding to this node
                c_vector<double,DIM> image_node_location = real_node_location;
                if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
                {
                    image_node_location[0] -= mPeriodicDomainWidth;
                }
                else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
                {
                    image_node_location[0] += mPeriodicDomainWidth;
                }

                // Create a copy of the node corresponding to this cell, suitable translated, and store it
                Node<DIM>* p_image_node = new Node<DIM>(num_cells+count, image_node_location);
                extended_nodes[num_cells+count] = p_image_node;

                // Populate mExtendedMeshNodeIndexMap
                mExtendedMeshNodeIndexMap[count] = real_node_index;
                mExtendedMeshNodeIndexMap[num_cells+count] = real_node_index;

                count++;
            }

            // Second, extend this extended mesh in the y-direction, so that we cover the corners too
            // (We don't need to store the real nodes anymore

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
            {
                // First, create and store a copy of this real node and cell
                unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                // Compute the location of the image node corresponding to this node
                c_vector<double,DIM> image_node_location = real_node_location;

                if (real_node_location[1] >= mPeriodicDomainDepth*0.5)
                {
                    image_node_location[1] -= mPeriodicDomainDepth;
                }
                else if (real_node_location[1] <  mPeriodicDomainDepth*0.5)
                {
                    image_node_location[1] += mPeriodicDomainDepth;
                }

                // Create a copy of the node corresponding to this cell, suitable translated, and store it
                Node<DIM>* p_image_node = new Node<DIM>(num_cells+count, image_node_location);
                extended_nodes[num_cells+count] = p_image_node;

                // Populate mExtendedMeshNodeIndexMap
                mExtendedMeshNodeIndexMap[num_cells+count] = real_node_index;

                count++;
            }

            // We now construct mpExtendedMesh using extended_nodes
            mpExtendedMesh = new MutableMesh<DIM,DIM>(extended_nodes);

        	// Now loop over the extended mesh and calculate the force acting on real nodes
        	// (using the edge iterator ensures that each edge is visited exactly once)
            for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator = mpExtendedMesh->EdgesBegin();
                 edge_iterator != mpExtendedMesh->EdgesEnd();
                 ++edge_iterator)
            {
                unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
                unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

                c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation);

                // Apply this force to any real nodes (i.e. nodes whose indices are less than num_real_nodes)
                if (nodeA_global_index < num_cells)
                {
                    unsigned real_node_index_A = mExtendedMeshNodeIndexMap[nodeA_global_index];
                    rForces[real_node_index_A] += force;
                }
                if (nodeB_global_index < num_cells)
                {
                    unsigned real_node_index_B = mExtendedMeshNodeIndexMap[nodeB_global_index];
                    rForces[real_node_index_B] -= force;
                }
            }

			break;
        }

        default:
            // This can't happen
            NEVER_REACHED;
    }

}

template<unsigned DIM>
double AbstractPeriodicTwoBodyInteractionForce<DIM>::GetPeriodicDomainWidth()
{
	return mPeriodicDomainWidth;
}

template<unsigned DIM>
void AbstractPeriodicTwoBodyInteractionForce<DIM>::SetPeriodicDomainWidth(double periodicDomainWidth)
{
	mPeriodicDomainWidth = periodicDomainWidth;
}

template<unsigned DIM>
double AbstractPeriodicTwoBodyInteractionForce<DIM>::GetPeriodicDomainDepth()
{
	return mPeriodicDomainDepth;
}

template<unsigned DIM>
void AbstractPeriodicTwoBodyInteractionForce<DIM>::SetPeriodicDomainDepth(double periodicDomainDepth)
{
	mPeriodicDomainDepth = periodicDomainDepth;
}

template<unsigned DIM>
void AbstractPeriodicTwoBodyInteractionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PeriodicDomainWidth>" << mPeriodicDomainWidth << "</PeriodicDomainWidth> \n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AbstractPeriodicTwoBodyInteractionForce<1>;
template class AbstractPeriodicTwoBodyInteractionForce<2>;
template class AbstractPeriodicTwoBodyInteractionForce<3>;
