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

template<unsigned DIM>
AbstractPeriodicTwoBodyInteractionForce<DIM>::AbstractPeriodicTwoBodyInteractionForce()
   : AbstractTwoBodyInteractionForce<DIM>(),
     mInitialWidth(0.0)
{
}

template<unsigned DIM>
void AbstractPeriodicTwoBodyInteractionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                                AbstractCellPopulation<DIM>& rCellPopulation)
{
    /*
     * Rather than iterating over the springs here, we first create a new, extended mesh by copying
     * real nodes to form image nodes on either side of the original mesh and creating an extended
     * cell population. We then use this extended cell population to work out the forces acting on
     * the real nodes.
     */

    // This method currently works only in 2d
    assert(DIM == 2);

    unsigned num_real_nodes = rCellPopulation.GetNumRealCells();
	std::vector<Node<DIM>*> extended_nodes(2*num_real_nodes);
    std::vector<CellPtr> extended_cells(2*num_real_nodes);

	// The width of the extended mesh
	double extended_mesh_width =  mInitialWidth;

	// Calculate the cell population's centroid
	c_vector<double, DIM> centroid = rCellPopulation.GetCentroidOfCellPopulation();

    ///\todo The code block below would need to be amended for 3d (#1856)
    ///\todo Think more about the use of centroid to make more general (#1856)

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        // Create a copy of this cell and store it
        CellPtr p_real_cell = *cell_iter;
        extended_cells[real_node_index] = p_real_cell;

        // Create a copy of the node corresponding to this cell and store it
    	Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
    	extended_nodes[real_node_index] = p_real_node;

        // Second, create and store the corresponding image node and cell
        unsigned image_node_index = real_node_index + num_real_nodes;
        c_vector<double, DIM> image_node_location = real_node_location;
        if (image_node_location[0] >= centroid(0))
        {
            image_node_location[0] -= extended_mesh_width;
        }
        else if (image_node_location[0] < centroid(0))
        {
            image_node_location[0] += extended_mesh_width;
        }

        // Create a copy of this cell and store it
        CellPtr p_image_cell = *cell_iter;
        extended_cells[image_node_index] = p_image_cell;

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<DIM>* p_image_node = new Node<DIM>(image_node_index, image_node_location);
        extended_nodes[image_node_index] = p_image_node;
    }

    // Now construct a mesh using extended_nodes
    MutableMesh<DIM,DIM>* p_extended_mesh = new MutableMesh<DIM,DIM>(extended_nodes);

    // Use extended_mesh and extended_cells to create a MeshBasedCellPopulation
    MeshBasedCellPopulation<DIM>* p_extended_cell_population = new MeshBasedCellPopulation<DIM>(*p_extended_mesh, extended_cells, std::vector<unsigned>(), false, false);

    ///\todo should we call p_extended_cell_population->Update() to ensure mMarkedSprings is correct? (#1856)

	// Now loop over the extended mesh and calculate the force acting on real nodes
	// (using the edge iterator ensures that each edge is visited exactly once)
    for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator = p_extended_mesh->EdgesBegin();
         edge_iterator != p_extended_mesh->EdgesEnd();
         ++edge_iterator)
    {
        unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, *p_extended_cell_population);

        // Apply this force to any real nodes (i.e. nodes whose indices are less than num_real_nodes)
        if (nodeA_global_index < num_real_nodes)
        {
            rForces[nodeA_global_index] += force;
        }
        if (nodeB_global_index < num_real_nodes)
        {
            rForces[nodeB_global_index] -= force;
        }
    }

    // Avoid memory leaks
    delete p_extended_cell_population;
    delete p_extended_mesh;
}

template<unsigned DIM>
double AbstractPeriodicTwoBodyInteractionForce<DIM>::GetInitialWidth()
{
	return mInitialWidth;
}

template<unsigned DIM>
void AbstractPeriodicTwoBodyInteractionForce<DIM>::SetInitialWidth(double initialWidth)
{
	mInitialWidth = initialWidth;
}

template<unsigned DIM>
void AbstractPeriodicTwoBodyInteractionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<InitialWidth>" << mInitialWidth << "</InitialWidth> \n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AbstractPeriodicTwoBodyInteractionForce<1>;
template class AbstractPeriodicTwoBodyInteractionForce<2>;
template class AbstractPeriodicTwoBodyInteractionForce<3>;
