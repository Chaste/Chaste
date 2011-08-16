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

#include "PottsMeshGenerator.hpp"
template<unsigned DIM>
PottsMeshGenerator<DIM>::PottsMeshGenerator( unsigned numNodesAcross, unsigned numElementsAcross, unsigned elementWidth,
											 unsigned numNodesUp, unsigned numElementsUp, unsigned elementHeight,
											 unsigned numNodesDeep, unsigned numElementsDeep, unsigned elementDepth,
											 bool startAtBottomLeft)
{
    assert(numElementsAcross > 0);
    assert(numElementsUp > 0);
    assert(elementWidth > 0);
    assert(elementHeight > 0);
    assert(numElementsDeep>0);
    assert(numNodesDeep > 0);
    assert(elementDepth>0);
    assert(numElementsAcross*elementWidth<=numNodesAcross);
    assert(numElementsUp*elementHeight<=numNodesUp);
    assert(numElementsDeep*elementDepth<=numNodesDeep);

    std::vector<Node<DIM>*> nodes;
    std::vector<PottsElement<DIM>*>  elements;

    unsigned node_index = 0;
    unsigned node_indices[elementWidth*elementHeight*elementDepth];
    unsigned element_index;

    // Calculate the width of the medium on the edge and offset the node index so that the elements are in the centre of the mesh.
    unsigned across_gap = (numNodesAcross -  numElementsAcross*elementWidth)/2;
    unsigned up_gap = (numNodesUp -  numElementsUp*elementHeight)/2;
    unsigned deep_gap = (numNodesDeep -  numElementsDeep*elementDepth)/2;

    unsigned index_offset = 0;

    if (!startAtBottomLeft) // Elements in centre of mesh
    {
        index_offset = deep_gap*numNodesAcross*numNodesUp +
                            up_gap*numNodesAcross +
                            across_gap;
    }

    /*
     * Create the nodes, row by row, from the bottom up
     * On the first and last row we have numNodesAcross nodes, all of which are boundary
     * nodes. On each interior row we have numNodesAcross nodes, the first and last nodes
     * are boundary nodes.
     */
    for (unsigned k=0; k<numNodesDeep; k++)
    {
		for (unsigned j=0; j<numNodesUp; j++)
		{
			for (unsigned i=0; i<numNodesAcross; i++)
			{
				bool is_boundary_node=false;
				if (DIM==2)
				{
					is_boundary_node = (j==0 || j==numNodesUp-1 || i==0 || i==numNodesAcross-1) ? true : false;
				}
				if (DIM==3)
				{
					is_boundary_node = (j==0 || j==numNodesUp-1 || i==0 || i==numNodesAcross-1 || k==0 || k==numNodesDeep-1) ? true : false;
				}
				Node<DIM>* p_node = new Node<DIM>(node_index, is_boundary_node, i, j, k);
				nodes.push_back(p_node);
				node_index++;
			}
		}
    }

    /*
     * Create the elements. The array node_indices contains the
     * global node indices, in increasing order.
     */
    for (unsigned n=0; n<numElementsDeep; n++)
    {
		for (unsigned j=0; j<numElementsUp; j++)
		{
			for (unsigned i=0; i<numElementsAcross; i++)
			{
				for (unsigned m=0; m<elementDepth; m++)
				{
					for (unsigned l=0; l<elementHeight; l++)
					{
						for (unsigned k=0; k<elementWidth; k++)
						{
							node_indices[m*elementHeight*elementWidth + l*elementWidth + k] = n*elementDepth*numNodesUp*numNodesAcross +
							                                                                  j*elementHeight*numNodesAcross +
							                                                                  i*elementWidth +
							                                                                  m*numNodesAcross*numNodesUp +
							                                                                  l*numNodesAcross +
							                                                                  k + index_offset;

						}
					}
				}
				std::vector<Node<DIM>*> element_nodes;
				for (unsigned k=0; k<elementDepth*elementHeight*elementWidth; k++)
				{
				   element_nodes.push_back(nodes[node_indices[k]]);
				}

				element_index = n*numElementsAcross*numElementsUp + j*numElementsAcross + i;
				PottsElement<DIM>* p_element = new PottsElement<DIM>(element_index, element_nodes);
				elements.push_back(p_element);
			}
		}
    }

    mpMesh = new PottsMesh<DIM>(nodes, elements);
}
template<unsigned DIM>
PottsMeshGenerator<DIM>::~PottsMeshGenerator()
{
    delete mpMesh;
}

template<unsigned DIM>
PottsMesh<DIM>* PottsMeshGenerator<DIM>::GetMesh()
{
    return mpMesh;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class PottsMeshGenerator<1>;
template class PottsMeshGenerator<2>;
template class PottsMeshGenerator<3>;
