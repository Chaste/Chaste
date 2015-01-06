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
#include "BoxCollection.hpp"
#include "Exception.hpp"

/////////////////////////////////////////////////////////////////////////////
// BoxCollection methods
/////////////////////////////////////////////////////////////////////////////

// Static member for "fudge factor" is instantiated here
template<unsigned DIM>
const double BoxCollection<DIM>::msFudge = 5e-14;

template<unsigned DIM>
BoxCollection<DIM>::BoxCollection(double boxWidth, c_vector<double, 2*DIM> domainSize, bool isPeriodicInX)
    : mDomainSize(domainSize),
      mBoxWidth(boxWidth),
      mIsPeriodicInX(isPeriodicInX)
{
    // Periodicity only works in 2d
    if (isPeriodicInX)
    {
        assert(DIM==2);
    }

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
        mNumBoxesEachDirection(i) = (unsigned) floor((domainSize(2*i+1) - domainSize(2*i))/boxWidth + msFudge) + 1;
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
void BoxCollection<DIM>::EmptyBoxes()
{
    for (unsigned i=0; i<mBoxes.size(); i++)
    {
        mBoxes[i].ClearNodes();
    }
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
        containing_box_indices[i] = (unsigned) floor((rLocation[i] - mDomainSize(2*i) + msFudge)/mBoxWidth);
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
                    // If we're on the left edge but its periodic include the box on the far right
                    else if ( (box_index % mNumBoxesEachDirection(0) == 0) && (mIsPeriodicInX) )
                    {
                        local_boxes.insert(box_index + 2* mNumBoxesEachDirection(0) - 1);
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
                // If we're on the right edge but it's periodic include the box on the far left of the domain
                else if ( (box_index % mNumBoxesEachDirection(0) == mNumBoxesEachDirection(0)-1) && (mIsPeriodicInX) )
                {
                    local_boxes.insert(box_index - mNumBoxesEachDirection(0) + 1);
                    // If we're also not on the top-most row, then insert the box above- on the far left of the domain
                    if (box_index < mBoxes.size() - mNumBoxesEachDirection(0))
                    {
                        local_boxes.insert(box_index + 1);
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
                // If we're not on the bottom face (z min), then DON'T insert the box above - this will lead to duplicate pairs.
                if (box_index >= num_boxes_xy)
                {
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
                else // Add Periodic Box if needed
                {
                    if(mIsPeriodicInX)
                    {
                        local_boxes.insert(i+M-1);
                    }
                }

                // add the box to the right
                if (!is_xmax[i])
                {
                    local_boxes.insert(i+1);
                }
                else // Add Periodic Box if needed
                {
                    if(mIsPeriodicInX)
                    {
                        local_boxes.insert(i-M+1);
                    }
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

                // Add Periodic Corner Boxes if needed
                if(mIsPeriodicInX)
                {
                    if( (is_xmin[i]) && (!is_ymin[i]) )
                    {
                        local_boxes.insert(i-1);
                    }
                    if ( (is_xmin[i]) && (!is_ymax[i]) )
                    {
                        local_boxes.insert(i-1+2*M);
                    }
                    if ( (is_xmax[i]) && (!is_ymin[i]) )
                    {
                        local_boxes.insert(i+1-2*M);
                    }
                    if ( (is_xmax[i]) && (!is_ymax[i]) )
                    {
                        local_boxes.insert(i+1);
                    }
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
const c_vector<double, 2*DIM>& BoxCollection<DIM>::rGetDomainSize() const
{
    return mDomainSize;
}

template<unsigned DIM>
void BoxCollection<DIM>::CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs, std::map<unsigned, std::set<unsigned> >& rNodeNeighbours)
{
    rNodePairs.clear();
    rNodeNeighbours.clear();

    // Create an empty set of neighbours for each node
    for (unsigned node_index=0; node_index<rNodes.size(); node_index++)
    {
        rNodeNeighbours[node_index] = std::set<unsigned>();
    }

    for (unsigned i=0; i<rNodes.size(); i++)
    {
        unsigned node_index = rNodes[i]->GetIndex();

        // Get the box containing this node
        unsigned box_index = CalculateContainingBox(rNodes[i]);

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
                    if (other_node_index > rNodes[i]->GetIndex())
                    {
                        rNodePairs.push_back(std::pair<Node<DIM>*, Node<DIM>*>(rNodes[i], (*node_iter)));
                        rNodeNeighbours[node_index].insert(other_node_index);
                        rNodeNeighbours[other_node_index].insert(node_index);
                    }
                }
                else
                {
                    rNodePairs.push_back(std::pair<Node<DIM>*, Node<DIM>*>(rNodes[i], (*node_iter)));
                    rNodeNeighbours[node_index].insert(other_node_index);
                    rNodeNeighbours[other_node_index].insert(node_index);
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class BoxCollection<1>;
template class BoxCollection<2>;
template class BoxCollection<3>;
