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
    // This method currently works only in 2d
    assert(DIM == 2);

    // Create a helper pointer
	MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

	// Rather than iterating over the springs here, we firstly want to extend the mesh by creating
	// the image nodes, then loop over all edges to work out the forces acting on the real nodes.

	// Create a new, extended mesh by copying real nodes to form image nodes on either side, and create an accompanying extended cell population

	// These vectors will contain the real nodes and image nodes, and real cells and image cells, respectively
	std::vector<Node<DIM>*> extended_node_set;
    std::vector<CellPtr > extended_cell_set;

	// This vector will contain only the image nodes, and image cells, respectively
	std::vector<Node<DIM>*> image_node_set;
    std::vector<CellPtr> image_cells;

    // This vector will contain only the real cells
	std::vector<CellPtr> real_cells;

	unsigned num_real_nodes = rCellPopulation.GetNumRealCells();
	unsigned new_image_node_index = num_real_nodes;

	// The width of the extended mesh
	double extended_mesh_width =  mInitialWidth;

	// Calculate centroid of cell population
	c_vector<double, DIM> centroid = rCellPopulation.GetCentroidOfCellPopulation();

    // We iterate over all cells in the population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Create the 'real' cell corresponding to this cell and add it to the vector of real cells
        CellPtr p_real_cell(new Cell(cell_iter->GetMutationState(), cell_iter->GetCellCycleModel(), false));
        real_cells.push_back(p_real_cell);

        // Get the node corresponding to this cell
    	Node<DIM>* p_node = p_cell_population->GetNodeCorrespondingToCell(*cell_iter);
    	c_vector<double,DIM> real_node_location = p_node->rGetLocation();
    	unsigned real_node_index = p_node->GetIndex();

    	// Create a copy of this node and add it to the vector of all nodes
    	// (we push back all the original nodes first so that they are all kept together)
    	Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
    	extended_node_set.push_back(p_real_node);

    	// todo The code block below would need to be amended for 3d
    	// todo Think more about the use of centroid to make more general

    	// Compute the location of the image node corresponding to this node
        c_vector<double,DIM> image_node_location = real_node_location;
		if (real_node_location[0] >= centroid(0)) // Right-hand boundary node
		{
			image_node_location[0] -= extended_mesh_width;
		}
		else if (real_node_location[0] < centroid(0))
		{
			image_node_location[0] += extended_mesh_width;
		}

		// Create the image node corresponding to this node and add it to the vector of image nodes
		Node<DIM>* p_image_node = new Node<DIM>(new_image_node_index, image_node_location);
        image_node_set.push_back(p_image_node);

		// Create the image cell corresponding to this cell and add it to the vector of image cells
		CellPtr p_image_cell(new Cell(cell_iter->GetMutationState(), cell_iter->GetCellCycleModel(), false));
		image_cells.push_back(p_image_cell);

		// todo No set age method, which may matter when it comes to working out the rest length
		//      of the spring for the force calculation
//		p_new_image_cell->SetAge();

        // Start from total number of real nodes and increment upwards
		new_image_node_index++;
    }

    // Now construct the vectors extended_node_set and extended_cell_set so that
    // the image nodes/cells are together at the end of each vector
    for (unsigned i=0; i<image_node_set.size(); i++)
    {
    	extended_node_set.push_back(image_node_set[i]);
    }
    for (unsigned i=0; i<real_cells.size(); i++)
    {
    	extended_cell_set.push_back(real_cells[i]);
    }
    for (unsigned i=0; i<image_cells.size(); i++)
    {
    	extended_cell_set.push_back(image_cells[i]);
    }

    // Check that the vectors extended_node_set and extended_cell_set are the correct size
    assert(extended_cell_set.size() == 2*num_real_nodes);
    assert(extended_node_set.size() == 2*num_real_nodes);

    // We now construct a mesh using extended_node_set...
    MutableMesh<DIM,DIM> extended_mesh(extended_node_set);

    // ...and, with this mesh and extended_cell_set, we create a MeshBasedCellPopulation
    MeshBasedCellPopulation<DIM>* p_extended_cell_population = new MeshBasedCellPopulation<DIM>(extended_mesh, extended_cell_set);

    // todo might we need to call Update() on extended_cell_population to ensure
    // that mMarkedSprings is correct?

	// Now loop over the extended mesh and calculate the force acting on real nodes
	// (using the edge iterator ensures that each edge is visited only once)
    for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator = extended_mesh.EdgesBegin();
         edge_iterator != extended_mesh.EdgesEnd();
         ++edge_iterator)
    {
        unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index =  edge_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, *p_extended_cell_population);

        // Now we make sure that we only apply the force to the real node and not the image node
        if ( (nodeA_global_index < num_real_nodes) &&  (nodeB_global_index < num_real_nodes) )
        {
			rForces[nodeB_global_index] -= force;
			rForces[nodeA_global_index] += force;
        }
        else if ( (nodeA_global_index >= num_real_nodes) &&  (nodeB_global_index < num_real_nodes) )
        {
			rForces[nodeB_global_index] -= force;
        }
        else if ( (nodeA_global_index < num_real_nodes) &&  (nodeB_global_index >= num_real_nodes) )
        {
        	rForces[nodeA_global_index] += force;
        }
    }
}

template<unsigned DIM>
double AbstractPeriodicTwoBodyInteractionForce<DIM>::GetInitialWidth()
{
	return mInitialWidth;
}

template<unsigned DIM>
void AbstractPeriodicTwoBodyInteractionForce<DIM>::SetInitialWidth(double initialWidth)
{
	mInitialWidth=initialWidth;
}

template<unsigned DIM>
void AbstractPeriodicTwoBodyInteractionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AbstractPeriodicTwoBodyInteractionForce<1>;
template class AbstractPeriodicTwoBodyInteractionForce<2>;
template class AbstractPeriodicTwoBodyInteractionForce<3>;
