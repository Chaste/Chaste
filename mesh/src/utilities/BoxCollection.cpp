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
#include "BoxCollection.hpp"
#include "Exception.hpp"

/////////////////////////////////////////////////////////////////////////////
// Box methods
/////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
Box<DIM>::Box(c_vector<double, 2*DIM>& rMinAndMaxValues)
{
    mMinAndMaxValues = rMinAndMaxValues;
}

template<unsigned DIM>
c_vector<double, 2*DIM>& Box<DIM>::rGetMinAndMaxValues()
{
    return mMinAndMaxValues;
}

template<unsigned DIM>
void Box<DIM>::AddNode(Node<DIM>* pNode)
{
    mNodesContained.insert(pNode);
}

template<unsigned DIM>
void Box<DIM>::RemoveNode(Node<DIM>* pNode)
{
    mNodesContained.erase(pNode);
}

template<unsigned DIM>
std::set< Node<DIM>* >& Box<DIM>::rGetNodesContained()
{
    return mNodesContained;
}

template<unsigned DIM>
void Box<DIM>::AddElement(Element<DIM,DIM>* pElement)
{
    mElementsContained.insert(pElement);
}

template<unsigned DIM>
std::set< Element<DIM,DIM>* >& Box<DIM>::rGetElementsContained()
{
    return mElementsContained;
}


/////////////////////////////////////////////////////////////////////////////
// BoxCollection methods
/////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
BoxCollection<DIM>::BoxCollection(double boxWidth, c_vector<double, 2*DIM> domainSize)
    : mDomainSize(domainSize),
      mBoxWidth(boxWidth)
{
    /*
     * Start by calculating the number of boxes in each direction and total number of boxes.
     * Also create a helper vector of coefficients, whose first entry is 1 and whose i-th
     * entry (for i>1) is the i-th partial product of the vector mNumBoxesEachDirection. This
     * vector of coefficients will be used in the next code block to compute how many boxes
     * along in each dimension each box, given its index.
     */
    unsigned num_boxes = 1;
    std::vector<unsigned> coefficients;
    coefficients.push_back(1);

    for (unsigned i=0; i<DIM; i++)
    {
        mNumBoxesEachDirection(i) = (unsigned)((domainSize(2*i+1) - domainSize(2*i))/boxWidth) + 1;
        num_boxes *= mNumBoxesEachDirection(i);
        coefficients.push_back(coefficients[i]*mNumBoxesEachDirection(i));
    }

    for (unsigned box_index=0; box_index<num_boxes; box_index++)
    {
        /*
         * The code block below computes how many boxes along in each dimension the
         * current box is and stores this information in the second, ..., (DIM+1)th
         * entries of the vector current_box_indices. The first entry of
         * current_box_indices is zero to ensure that the for loop works correctly.
         *
         * Our convention is that in 3D the index of each box, box_index, is related
         * to its indices (i,j,k) by
         *
         * box_index = i + mNumBoxesEachDirection(0)*j
         *               + mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1)*k,
         *
         * while in 2D, box_index is related to (i,j) by
         *
         * box_index = i + mNumBoxesEachDirection(0)*j
         *
         * and in 1D we simply have box_index = i.
         */
        c_vector<unsigned, DIM+1> current_box_indices;
        current_box_indices[0] = 0;

        for (unsigned i=0; i<DIM; i++)
        {
            unsigned temp = 0;
            for (unsigned j=1; j<i; j++)
            {
                temp += coefficients[j]*current_box_indices[j-1];
            }
            current_box_indices[i+1] = (box_index%coefficients[i+1] - temp)/coefficients[i];
        }

        /*
         * We now use the information stores in current_box_indices to construct the
         * Box, which we add to mBoxes.
         */
        c_vector<double, 2*DIM> box_coords;
        for (unsigned l=0; l<DIM; l++)
        {
            box_coords(2*l) = domainSize(2*l) + current_box_indices(l+1)*boxWidth;
            box_coords(2*l+1) = domainSize(2*l) + (current_box_indices(l+1)+1)*boxWidth;
        }

        Box<DIM> new_box(box_coords);
        mBoxes.push_back(new_box);
    }

    // Check that we have the correct number of boxes
    assert(num_boxes == mBoxes.size());
}

template<unsigned DIM>
unsigned BoxCollection<DIM>::CalculateContainingBox(Node<DIM>* pNode)
{
    // Get the location of the node
    c_vector<double, DIM> location = pNode->rGetLocation();
    return CalculateContainingBox(location);
}


template<unsigned DIM>
unsigned BoxCollection<DIM>::CalculateContainingBox(c_vector<double, DIM>& rLocation)
{
    // The node must lie inside the boundary of the box collection
    for (unsigned i=0; i<DIM; i++)
    {
        if ( (rLocation[i] < mDomainSize(2*i)) || (rLocation[i] > mDomainSize(2*i+1)) )
        {
            EXCEPTION("The point provided is outside all of the boxes");
        }
    }

    // Compute the containing box index in each dimension
    c_vector<unsigned, DIM> containing_box_indices;
    for (unsigned i=0; i<DIM; i++)
    {
        containing_box_indices[i] = (unsigned) floor((rLocation[i] - mDomainSize(2*i))/mBoxWidth);
    }

    // Use these to compute the index of the containing box
    unsigned containing_box_index = 0;
    for (unsigned i=0; i<DIM; i++)
    {
        unsigned temp = 1;
        for (unsigned j=0; j<i; j++)
        {
            temp *= mNumBoxesEachDirection(j);
        }
        containing_box_index += temp*containing_box_indices[i];
    }

    // This index must be less than the number of boxes
    assert(containing_box_index < mBoxes.size());

    return containing_box_index;
}

template<unsigned DIM>
Box<DIM>& BoxCollection<DIM>::rGetBox(unsigned boxIndex)
{
    assert(boxIndex < mBoxes.size());
    return mBoxes[boxIndex];
}

template<unsigned DIM>
unsigned BoxCollection<DIM>::GetNumBoxes()
{
    return mBoxes.size();
}

template<unsigned DIM>
void BoxCollection<DIM>::SetupLocalBoxesHalfOnly()
{
    switch (DIM)
    {
        case 1:
        {
            // We only need to look for neighbours in the current and successive boxes
            mLocalBoxes.clear();
            for (unsigned box_index=0; box_index<mBoxes.size(); box_index++)
            {
                std::set<unsigned> local_boxes;

                // Insert the current box
                local_boxes.insert(box_index);

                // If we're not at the right-most box, then insert the box to the right
                if (box_index < mNumBoxesEachDirection(0)-1)
                {
                    local_boxes.insert(box_index+1);
                }
                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        case 2:
        {
            // We only need to look for neighbours in the current box and half the neighbouring boxes
            mLocalBoxes.clear();
            for (unsigned box_index=0; box_index<mBoxes.size(); box_index++)
            {
                std::set<unsigned> local_boxes;

                // Insert the current box
                local_boxes.insert(box_index);

                // If we're not on the top-most row, then insert the box above
                if (box_index < mBoxes.size() - mNumBoxesEachDirection(0))
                {
                    local_boxes.insert(box_index + mNumBoxesEachDirection(0));

                    // If we're also not on the left-most column, then insert the box above-left
                    if (box_index % mNumBoxesEachDirection(0) != 0)
                    {
                        local_boxes.insert(box_index + mNumBoxesEachDirection(0) - 1);
                    }
                }
                // If we're not on the right-most column, then insert the box to the right
                if (box_index % mNumBoxesEachDirection(0) != mNumBoxesEachDirection(0)-1)
                {
                    local_boxes.insert(box_index + 1);

                    // If we're also not on the top-most row, then insert the box above-right
                    if (box_index < mBoxes.size() - mNumBoxesEachDirection(0))
                    {
                        local_boxes.insert(box_index + mNumBoxesEachDirection(0) + 1);
                    }
                }

                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        case 3:
        {
            // We only need to look for neighbours in the current box and half the neighbouring boxes
            mLocalBoxes.clear();
            unsigned num_boxes_xy = mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1);

            for (unsigned box_index=0; box_index<mBoxes.size(); box_index++)
            {
                std::set<unsigned> local_boxes;

                // Insert the current box
                local_boxes.insert(box_index);

                // If we're not on the far face (y max), then insert the far box
                if (box_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
                {
                    local_boxes.insert(box_index + mNumBoxesEachDirection(0));

                    // If we're also not on the left face (x min), then insert the box to the left
                    if (box_index % mNumBoxesEachDirection(0) != 0)
                    {
                        local_boxes.insert(box_index + mNumBoxesEachDirection(0) - 1);
                    }
                }
                // If we're not on the right face (x max), then insert the box to the right
                if (box_index % mNumBoxesEachDirection(0) != mNumBoxesEachDirection(0)-1)
                {
                    local_boxes.insert(box_index + 1);

                    // If we're also not on the far face (y max) row, then insert the box to the far-right
                    if (box_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
                    {
                        local_boxes.insert(box_index + mNumBoxesEachDirection(0) + 1);
                    }
                }
                // If we're not on the top face (z max), then insert the box above
                if (box_index < mBoxes.size() - num_boxes_xy)
                {
                    local_boxes.insert(box_index + num_boxes_xy);

                    // If we're also not on the far face (y max), then insert the above-far box
                    if (box_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
                    {
                        local_boxes.insert(box_index + num_boxes_xy + mNumBoxesEachDirection(0));

                        // If we're also not on the left face (x min), then insert the box to the above-left
                        if (box_index % mNumBoxesEachDirection(0) != 0)
                        {
                            local_boxes.insert(box_index + num_boxes_xy + mNumBoxesEachDirection(0) - 1);
                        }
                    }
                    // If we're also not on the right face (x max), then insert the box to the above-right
                    if (box_index % mNumBoxesEachDirection(0) != mNumBoxesEachDirection(0)-1)
                    {
                        local_boxes.insert(box_index + num_boxes_xy + 1);

                        // If we're also not on the far face (y max) row, then insert the box to the above-far-right
                        if (box_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
                        {
                            local_boxes.insert(box_index + num_boxes_xy + mNumBoxesEachDirection(0) + 1);
                        }
                    }
                }
                // If we're not on the bottom face (z min), then insert the box above
                if (box_index >= num_boxes_xy)
                {
                    local_boxes.insert(box_index - num_boxes_xy);

                    // If we're also not on the far face (y max), then insert the below-far box
                    if (box_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
                    {
                        local_boxes.insert(box_index - num_boxes_xy + mNumBoxesEachDirection(0));

                        // If we're also not on the left face (x min), then insert the box to the below-left
                        if (box_index % mNumBoxesEachDirection(0) != 0)
                        {
                            local_boxes.insert(box_index - num_boxes_xy + mNumBoxesEachDirection(0) - 1);
                        }
                    }
                    // If we're also not on the right face (x max), then insert the box to the below-right
                    if (box_index % mNumBoxesEachDirection(0) != mNumBoxesEachDirection(0)-1)
                    {
                        local_boxes.insert(box_index - num_boxes_xy + 1);

                        // If we're also not on the far face (y max) row, then insert the box to the below-far-right
                        if (box_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0))
                        {
                            local_boxes.insert(box_index - num_boxes_xy + mNumBoxesEachDirection(0) + 1);
                        }
                    }
                }
                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        default:
            NEVER_REACHED;
    }
}



template<unsigned DIM>
void BoxCollection<DIM>::SetupAllLocalBoxes()
{
    switch (DIM)
    {
        case 1:
        {
            for (unsigned i=0; i<mBoxes.size(); i++)
            {
                std::set<unsigned> local_boxes;

                local_boxes.insert(i);

                // add the two neighbours
                if (i!=0)
                {
                    local_boxes.insert(i-1);
                }
                if (i+1 != mNumBoxesEachDirection(0))
                {
                    local_boxes.insert(i+1);
                }

                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        case 2:
        {
            mLocalBoxes.clear();

            unsigned M = mNumBoxesEachDirection(0);
            unsigned N = mNumBoxesEachDirection(1);

            std::vector<bool> is_xmin(N*M); // far left
            std::vector<bool> is_xmax(N*M); // far right
            std::vector<bool> is_ymin(N*M); // bottom
            std::vector<bool> is_ymax(N*M); // top

            for (unsigned i=0; i<M*N; i++)
            {
                is_xmin[i] = (i%M==0);
                is_xmax[i] = ((i+1)%M==0);
                is_ymin[i] = (i%(M*N)<M);
                is_ymax[i] = (i%(M*N)>=(N-1)*M);
            }

            for (unsigned i=0; i<mBoxes.size(); i++)
            {
                std::set<unsigned> local_boxes;

                local_boxes.insert(i);

                // add the box to the left
                if (!is_xmin[i])
                {
                    local_boxes.insert(i-1);
                }

                // add the box to the right
                if (!is_xmax[i])
                {
                    local_boxes.insert(i+1);
                }

                // add the one below
                if (!is_ymin[i])
                {
                    local_boxes.insert(i-M);
                }

                // add the one above
                if (!is_ymax[i])
                {
                    local_boxes.insert(i+M);
                }

                // add the four corner boxes

                if ( (!is_xmin[i]) && (!is_ymin[i]) )
                {
                    local_boxes.insert(i-1-M);
                }

                if ( (!is_xmin[i]) && (!is_ymax[i]) )
                {
                    local_boxes.insert(i-1+M);
                }

                if ( (!is_xmax[i]) && (!is_ymin[i]) )
                {
                    local_boxes.insert(i+1-M);
                }

                if ( (!is_xmax[i]) && (!is_ymax[i]) )
                {
                    local_boxes.insert(i+1+M);
                }

                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        case 3:
        {
            mLocalBoxes.clear();

            unsigned M = mNumBoxesEachDirection(0);
            unsigned N = mNumBoxesEachDirection(1);
            unsigned P = mNumBoxesEachDirection(2);

            std::vector<bool> is_xmin(N*M*P); // far left
            std::vector<bool> is_xmax(N*M*P); // far right
            std::vector<bool> is_ymin(N*M*P); // nearest
            std::vector<bool> is_ymax(N*M*P); // furthest
            std::vector<bool> is_zmin(N*M*P); // bottom layer
            std::vector<bool> is_zmax(N*M*P); // top layer

            for (unsigned i=0; i<M*N*P; i++)
            {
                is_xmin[i] = (i%M==0);
                is_xmax[i] = ((i+1)%M==0);
                is_ymin[i] = (i%(M*N)<M);
                is_ymax[i] = (i%(M*N)>=(N-1)*M);
                is_zmin[i] = (i<M*N);
                is_zmax[i] = (i>=M*N*(P-1));
            }

            for (unsigned i=0; i<mBoxes.size(); i++)
            {
                std::set<unsigned> local_boxes;

                // add itself as a local box
                local_boxes.insert(i);

                // now add all 26 other neighbours.....

                // add the box left
                if (!is_xmin[i])
                {
                    local_boxes.insert(i-1);

                    // plus some others towards the left
                    if (!is_ymin[i])
                    {
                        local_boxes.insert(i-1-M);
                    }

                    if (!is_ymax[i])
                    {
                        local_boxes.insert(i-1+M);
                    }

                    if (!is_zmin[i])
                    {
                        local_boxes.insert(i-1-M*N);
                    }

                    if (!is_zmax[i])
                    {
                        local_boxes.insert(i-1+M*N);
                    }
                }

                // add the box to the right
                if (!is_xmax[i])
                {
                    local_boxes.insert(i+1);

                    // plus some others towards the right
                    if (!is_ymin[i])
                    {
                        local_boxes.insert(i+1-M);
                    }

                    if (!is_ymax[i])
                    {
                        local_boxes.insert(i+1+M);
                    }

                    if (!is_zmin[i])
                    {
                        local_boxes.insert(i+1-M*N);
                    }

                    if (!is_zmax[i])
                    {
                        local_boxes.insert(i+1+M*N);
                    }
                }

                // add the boxes next along the y axis
                if (!is_ymin[i])
                {
                    local_boxes.insert(i-M);

                    // and more in this plane
                    if (!is_zmin[i])
                    {
                        local_boxes.insert(i-M-M*N);
                    }

                    if (!is_zmax[i])
                    {
                        local_boxes.insert(i-M+M*N);
                    }
                }

                // add the boxes next along the y axis
                if (!is_ymax[i])
                {
                    local_boxes.insert(i+M);

                    // and more in this plane
                    if (!is_zmin[i])
                    {
                        local_boxes.insert(i+M-M*N);
                    }

                    if (!is_zmax[i])
                    {
                        local_boxes.insert(i+M+M*N);
                    }
                }

                // add the box directly above
                if (!is_zmin[i])
                {
                    local_boxes.insert(i-N*M);
                }

                // add the box directly below
                if (!is_zmax[i])
                {
                    local_boxes.insert(i+N*M);
                }

                // finally, the 8 corners are left

                if ( (!is_xmin[i]) && (!is_ymin[i]) && (!is_zmin[i]) )
                {
                    local_boxes.insert(i-1-M-M*N);
                }

                if ( (!is_xmin[i]) && (!is_ymin[i]) && (!is_zmax[i]) )
                {
                    local_boxes.insert(i-1-M+M*N);
                }

                if ( (!is_xmin[i]) && (!is_ymax[i]) && (!is_zmin[i]) )
                {
                    local_boxes.insert(i-1+M-M*N);
                }

                if ( (!is_xmin[i]) && (!is_ymax[i]) && (!is_zmax[i]) )
                {
                    local_boxes.insert(i-1+M+M*N);
                }

                if ( (!is_xmax[i]) && (!is_ymin[i]) && (!is_zmin[i]) )
                {
                    local_boxes.insert(i+1-M-M*N);
                }

                if ( (!is_xmax[i]) && (!is_ymin[i]) && (!is_zmax[i]) )
                {
                    local_boxes.insert(i+1-M+M*N);
                }

                if ( (!is_xmax[i]) && (!is_ymax[i]) && (!is_zmin[i]) )
                {
                    local_boxes.insert(i+1+M-M*N);
                }

                if ( (!is_xmax[i]) && (!is_ymax[i]) && (!is_zmax[i]) )
                {
                    local_boxes.insert(i+1+M+M*N);
                }

                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        default:
            NEVER_REACHED;
    }
}

template<unsigned DIM>
std::set<unsigned> BoxCollection<DIM>::GetLocalBoxes(unsigned boxIndex)
{
    assert(boxIndex < mLocalBoxes.size());
    return mLocalBoxes[boxIndex];
}

template<unsigned DIM>
void BoxCollection<DIM>::CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::set<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    rNodePairs.clear();
    for (unsigned node_index=0; node_index<rNodes.size(); node_index++)
    {
        // Get the box containing this node
        unsigned box_index = CalculateContainingBox(rNodes[node_index]);

        // Get the local boxes to this node
        std::set<unsigned> local_boxes_indices = GetLocalBoxes(box_index);

        // Loop over all the local boxes
        for (std::set<unsigned>::iterator box_iter = local_boxes_indices.begin();
             box_iter != local_boxes_indices.end();
             box_iter++)
        {
            // Get the set of nodes contained in this box
            std::set< Node<DIM>* >& r_contained_nodes = mBoxes[*box_iter].rGetNodesContained();

            // Loop over these nodes
            for (typename std::set<Node<DIM>*>::iterator node_iter = r_contained_nodes.begin();
                 node_iter != r_contained_nodes.end();
                 ++node_iter)
            {
                // Get the index of the other node
                unsigned other_node_index = (*node_iter)->GetIndex();

                // If we're in the same box, then take care not to store the node pair twice
                if (*box_iter == box_index)
                {
                    if (other_node_index > node_index)
                    {
                        rNodePairs.insert(std::pair<Node<DIM>*, Node<DIM>*>(rNodes[node_index], rNodes[other_node_index]));
                    }
                }
                else
                {
                    rNodePairs.insert(std::pair<Node<DIM>*, Node<DIM>*>(rNodes[node_index], rNodes[other_node_index]));
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class Box<1>;
template class Box<2>;
template class Box<3>;
template class BoxCollection<1>;
template class BoxCollection<2>;
template class BoxCollection<3>;
