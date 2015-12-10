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
#include "Debug.hpp"

/////////////////////////////////////////////////////////////////////////////
// BoxCollection methods
/////////////////////////////////////////////////////////////////////////////

// Static member for "fudge factor" is instantiated here
template<unsigned DIM>
const double BoxCollection<DIM>::msFudge = 5e-14;

template<unsigned DIM>
BoxCollection<DIM>::BoxCollection(double boxWidth, c_vector<double, 2 * DIM> domainSize, bool isPeriodicInX,
                                  bool isPeriodicInY, bool isPeriodicInZ)
        : mDomainSize(domainSize),
          mBoxWidth(boxWidth)
{
    // Populate mIsDomainPeriodic
    switch (DIM)
    {
        case 1:
        {
            mIsDomainPeriodic[0] = isPeriodicInX;
            break;
        }
        case 2:
        {
            mIsDomainPeriodic[0] = isPeriodicInX;
            mIsDomainPeriodic[1] = isPeriodicInY;
            break;
        }
        case 3:
        {
            mIsDomainPeriodic[0] = isPeriodicInX;
            mIsDomainPeriodic[1] = isPeriodicInY;
            mIsDomainPeriodic[2] = isPeriodicInZ;
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }

    /*
     * Start by calculating the number of boxes in each direction and total number of boxes.
     * Also create a helper vector of coefficients, whose first entry is 1 and whose i-th
     * entry (for i>1) is the i-th partial product of the vector mNumBoxesEachDirection. This
     * vector of coefficients will be used in the next code block to compute how many boxes
     * along in each dimension each box is, given its index.
     */
    unsigned num_boxes = 1;
    std::vector<unsigned> coefficients;
    coefficients.push_back(1);

    for (unsigned i = 0; i < DIM; i++)
    {
        ///\todo #2725 example: domain width of 1.0 and box width of 0.25, the following line will create 5 boxes not 4
        mNumBoxesEachDirection(i) = (unsigned) floor((domainSize(2 * i + 1) - domainSize(2 * i)) / boxWidth + msFudge) + 1;
        num_boxes *= mNumBoxesEachDirection(i);
        coefficients.push_back(coefficients[i] * mNumBoxesEachDirection(i));
    }

    for (unsigned box_index = 0 ; box_index < num_boxes ; box_index++)
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
        c_vector<unsigned, DIM + 1> current_box_indices;
        current_box_indices[0] = 0;

        for (unsigned i = 0; i < DIM; i++)
        {
            unsigned temp = 0;
            for (unsigned j = 1; j < i; j++)
            {
                temp += coefficients[j] * current_box_indices[j - 1];
            }
            current_box_indices[i + 1] = (box_index % coefficients[i + 1] - temp) / coefficients[i];
        }

        /*
         * We now use the information stores in current_box_indices to construct the
         * Box, which we add to mBoxes.
         */
        c_vector<double, 2 * DIM> box_coords;
        for (unsigned l = 0; l < DIM; l++)
        {
            box_coords(2 * l) = domainSize(2 * l) + current_box_indices(l + 1) * boxWidth;
            box_coords(2 * l + 1) = domainSize(2 * l) + (current_box_indices(l + 1) + 1) * boxWidth;
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
    for (unsigned i = 0; i < mBoxes.size(); i++)
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
    // Compute the containing box index in each dimension
    c_vector<int, DIM> containing_box_indices;
    for (unsigned i = 0; i < DIM; i++)
    {
        containing_box_indices[i] = (int) floor((rLocation[i] - mDomainSize(2 * i) + msFudge) / mBoxWidth);
    }

    if(!IsBoxInDomain(containing_box_indices))
    {
        EXCEPTION("Location does not correspond to any box.");
    }

    return GetLinearIndex(containing_box_indices);
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
unsigned BoxCollection<DIM>::GetLinearIndex(c_vector<int, DIM> gridIndices)
{
    /*
     * This function may be passed values outside the range in one or more
     * dimensions in the case of a periodic domain.  We therefore assume that
     * these values represent a situation with periodicity and adjust them
     * accordingly before calculating the linear index.
     */

    // Adjust for periodicity if necessary
    for (unsigned dim = 0; dim < DIM; dim++)
    {
        // Check for values too large
        if (gridIndices(dim) >= (int)mNumBoxesEachDirection(dim))
        {
            assert(mIsDomainPeriodic(dim));
            gridIndices(dim) -= (int)mNumBoxesEachDirection(dim);
        }

        // Check for negative values
        else if (gridIndices(dim) < 0)
        {
            assert(mIsDomainPeriodic(dim));
            gridIndices(dim) += (int)mNumBoxesEachDirection(dim);
        }
    }

    // Calculate linear index
    unsigned linear_index;

    switch (DIM)
    {
        case 1:
        {
            linear_index = gridIndices(0);
            break;
        }
        case 2:
        {
            linear_index = gridIndices(0) +
                           gridIndices(1) * mNumBoxesEachDirection(0);
            break;
        }
        case 3:
        {
            linear_index = gridIndices(0) +
                           gridIndices(1) * mNumBoxesEachDirection(0) +
                           gridIndices(2) * mNumBoxesEachDirection(0) * mNumBoxesEachDirection(1);
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }

    return linear_index;
}

template<unsigned DIM>
c_vector<int,DIM> BoxCollection<DIM>::GetGridIndices(unsigned linearIndex)
{
    c_vector<int,DIM> grid_indices;

    switch (DIM)
    {
        case 1:
        {
            grid_indices(0) = linearIndex;
            break;
        }
        case 2:
        {
            unsigned num_x = mNumBoxesEachDirection(0);
            grid_indices(0) = linearIndex % num_x;
            grid_indices(1) = (linearIndex - grid_indices(0)) / num_x;
            break;
        }
        case 3:
        {
            unsigned num_x = mNumBoxesEachDirection(0);
            unsigned num_xy = mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1);
            grid_indices(0) = linearIndex % num_x;
            grid_indices(1) = (linearIndex % num_xy - grid_indices(0)) / num_x;
            grid_indices(2) = linearIndex / num_xy;
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }

    return grid_indices;
}

template<unsigned DIM>
bool BoxCollection<DIM>::IsBoxInDomain(c_vector<int, DIM> gridIndices)
{
    /*
     * We assume that, for a given dimension, any location is in the domain
     * if the domain is periodic in that dimension.
     *
     * This method can only return false in the case of one or more
     * dimensions being aperiodic and a location in one of those dimensions
     * being outside the range.
     */

    bool is_in_domain = true;

    for (unsigned dim = 0; dim < DIM; dim++)
    {
        if (!mIsDomainPeriodic(dim))
        {
            if (gridIndices(dim) < 0 || gridIndices(dim) >= (int)mNumBoxesEachDirection(dim))
            {
                is_in_domain = false;
            }
        }
    }

    return is_in_domain;
}

template<unsigned DIM>
c_vector<bool,DIM> BoxCollection<DIM>::IsIndexPenultimate(c_vector<int,DIM> gridIndices)
{
    c_vector<bool,DIM> is_penultimate;

    for (unsigned dim = 0 ; dim < DIM ; dim++)
    {
        is_penultimate(dim) = (gridIndices(dim) == (int)mNumBoxesEachDirection(dim) - 2);
    }

    return is_penultimate;
}

template<unsigned DIM>
void BoxCollection<DIM>::SetupLocalBoxesHalfOnly()
{
    // Populate a list of half the neighbours in this number of dimensions
    std::vector<c_vector<int,DIM> > neighbours;

    switch (DIM)
    {
        case 1:
        {
            // Just one neighbour, plus the current box (zero vector)
            neighbours = std::vector<c_vector<int,DIM> >(2);

            neighbours[0](0) = 0; // current box
            neighbours[1](0) = 1; // right

            break;
        }
        case 2:
        {
            // Four neighbours, plus the current box (zero vector)
            neighbours = std::vector<c_vector<int,DIM> >(5);

            neighbours[0](0) = 0; neighbours[0](1) = 1; // up
            neighbours[1](0) = 1; neighbours[1](1) = 1; // up right
            neighbours[2](0) = 1; neighbours[2](1) = 0; // right
            neighbours[3](0) = 1; neighbours[3](1) =-1; // down right
            neighbours[4](0) = 0; neighbours[4](1) = 0; // current box

            break;
        }
        case 3:
        {
            /*
             * We need to pick 13 neighbours in such a way as to ensure all interactions are captured,
             * in addition to the current box.
             *
             * The 26 cubes on the outside of a 3x3 arrangement either adjacent to a whole face of
             * the central cube, only an edge of the central cube, or only a vertex of the central
             * cube:
             *
             *    6 contact faces, and come in the following 3 pairs:
             *        +x        -x                  (1, 0, 0)
             *        +y        -y                  (0, 1, 0)
             *        +z        -z                  (0, 0, 1)
             *
             *    12 contact edges, and come in the following 6 pairs:
             *        +x+y      -x-y                (1, 1, 0)
             *        +x+z      -x-z                (1, 0, 1)
             *        +x-y      -x+y                (1,-1, 0)
             *        +x-z      -x+z                (1, 0,-1)
             *        +y+z      -y-z                (0, 1, 1)
             *        +y-z      -y+z                (0, 1,-1)
             *
             *    8 contact vertices, and come in the following 4 pairs:
             *        +x+y+z    -x-y-z              (1, 1, 1)
             *        +x+y-z    -x-y+z              (1, 1,-1)
             *        +x-y+z    -x+y-z              (1,-1, 1)
             *        +x-y-z    -x+y+z              (1,-1,-1)
             *
             * We will simply pick the 13 from the left-hand column above, which can be represented
             * as an integer offset in three dimensions from the middle cube, as shown in the third
             * column above.
             */

            neighbours = std::vector<c_vector<int,DIM> >(14);

            neighbours[0](0)  = 1; neighbours[0](1)  = 0; neighbours[0](2)  = 0;
            neighbours[1](0)  = 0; neighbours[1](1)  = 1; neighbours[1](2)  = 0;
            neighbours[2](0)  = 0; neighbours[2](1)  = 0; neighbours[2](2)  = 1;
            neighbours[3](0)  = 1; neighbours[3](1)  = 1; neighbours[3](2)  = 0;
            neighbours[4](0)  = 1; neighbours[4](1)  = 0; neighbours[4](2)  = 1;
            neighbours[5](0)  = 1; neighbours[5](1)  =-1; neighbours[5](2)  = 0;
            neighbours[6](0)  = 1; neighbours[6](1)  = 0; neighbours[6](2)  =-1;
            neighbours[7](0)  = 0; neighbours[7](1)  = 1; neighbours[7](2)  = 1;
            neighbours[8](0)  = 0; neighbours[8](1)  = 1; neighbours[8](2)  =-1;
            neighbours[9](0)  = 1; neighbours[9](1)  = 1; neighbours[9](2)  = 1;
            neighbours[10](0) = 1; neighbours[10](1) = 1; neighbours[10](2) =-1;
            neighbours[11](0) = 1; neighbours[11](1) =-1; neighbours[11](2) = 1;
            neighbours[12](0) = 1; neighbours[12](1) =-1; neighbours[12](2) =-1;

            neighbours[13](0) = 0; neighbours[13](1) = 0; neighbours[13](2) = 0; // current box

            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }

    // Pass the list of possible neighbours to SetupLocalBoxes()
    SetupLocalBoxes(neighbours);
}

template<unsigned DIM>
void BoxCollection<DIM>::SetupAllLocalBoxes()
{
    // Populate a list of all neighbours in this number of dimensions
    std::vector<c_vector<int,DIM> > neighbours;

    switch (DIM)
    {
        case 1:
        {
            // 3 neighbours
            neighbours.clear();

            for (int x_dim = -1 ; x_dim < 2 ; x_dim++)
            {
                c_vector<int,DIM> neighbour;
                neighbour(0) = x_dim;

                neighbours.push_back(neighbour);
            }

            break;
        }
        case 2:
        {
            // 3x3 neighbours
            neighbours.clear();

            for (int x_dim = -1 ; x_dim < 2 ; x_dim++)
            {
                for (int y_dim = -1 ; y_dim < 2 ; y_dim++)
                {
                    c_vector<int,DIM> neighbour;
                    neighbour(0) = x_dim;
                    neighbour(1) = y_dim;

                    neighbours.push_back(neighbour);
                }
            }

            break;
        }
        case 3:
        {
            // 3x3x3 neighbours
            neighbours.clear();

            for (int x_dim = -1 ; x_dim < 2 ; x_dim++)
            {
                for (int y_dim = -1 ; y_dim < 2 ; y_dim++)
                {
                    for (int z_dim = -1 ; z_dim < 2 ; z_dim++)
                    {
                        c_vector<int,DIM> neighbour;
                        neighbour(0) = x_dim;
                        neighbour(1) = y_dim;
                        neighbour(2) = z_dim;

                        neighbours.push_back(neighbour);
                    }
                }
            }

            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }

    // Pass the list of possible neighbours to SetupLocalBoxes()
    SetupLocalBoxes(neighbours);
}

template<unsigned DIM>
void BoxCollection<DIM>::SetupLocalBoxes(std::vector<c_vector<int, DIM> > neighbours)
{
    mLocalBoxes.clear();

    // Loop over all boxes and add all necessary local boxes
    for (unsigned box_idx = 0 ; box_idx < mBoxes.size() ; box_idx++)
    {
        std::set<unsigned> local_boxes;

        // The grid indices (i), (i,j) or (i,j,k) depending on DIM, corresponding to box_idx
        c_vector<int, DIM> grid_indices = GetGridIndices(box_idx);

        // Check all neighbours (1 or 2 in 1D, 4 or 8 in 2D, 13 or 26 in 3D) and add them as local
        // boxes if they are in the domain
        for (unsigned neighbour = 0 ; neighbour < neighbours.size() ; neighbour++)
        {
            if (IsBoxInDomain(grid_indices + neighbours[neighbour]))
            {
                local_boxes.insert(GetLinearIndex(grid_indices + neighbours[neighbour]));
            }
        }

        /*
         * If we are in a periodic domain, we need to add additional boxes if the
         * current box is in the penultimate position of any dimension.
         *
         * The additional boxes we need to add are just all the neighbours of the box
         * in the final position in that dimension.  Because we use set::insert, it
         * doesn't matter if we try and add the same box multiple times (only one
         * copy of the index will be included).
         */

        // Find if the current indices are penultimate in any dimension in which collection is periodic
        c_vector<bool,DIM> is_penultimate_periodic = IsIndexPenultimate(grid_indices);
        unsigned num_penultimate_periodic_dimensions = 0;

        for (unsigned dim = 0 ; dim < DIM ; dim++)
        {
            is_penultimate_periodic(dim) = is_penultimate_periodic(dim) && mIsDomainPeriodic(dim);
            if (is_penultimate_periodic(dim))
            {
                num_penultimate_periodic_dimensions++;
            }
        }

        // First, deal with at least one penultimate dimension
        if (num_penultimate_periodic_dimensions > 0)
        {
            // Loop over each dimension.  If penultimate in dimension DIM move one box and add all neighbours
            for (unsigned dim = 0 ; dim < DIM ; dim++)
            {
                if (is_penultimate_periodic(dim))
                {
                    // If we're penultimate in dimension dim, move to the final location and add each neighbour as before.
                    // The call to IsBoxInDomain() will handle periodicity.
                    c_vector<int,DIM> ultimate_indices = grid_indices;
                    ultimate_indices(dim) ++;

                    for (unsigned neighbour = 0 ; neighbour < neighbours.size() ; neighbour++)
                    {
                        if (IsBoxInDomain(ultimate_indices + neighbours[neighbour]))
                        {
                            local_boxes.insert(GetLinearIndex(ultimate_indices + neighbours[neighbour]));
                        }
                    }
                }
            }
        }

        // The final consideration is cases of multiple dimensions of periodicity.
        if (num_penultimate_periodic_dimensions > 1)
        {
            switch (DIM)
            {
                case 2:
                {
                    // Must have x and y periodicity - just have to add one extra box
                    assert(mIsDomainPeriodic(0) && mIsDomainPeriodic(1));

                    local_boxes.insert(0);

                    break;
                }
                case 3:
                {
                    // Could have X and Y, X and Z, Y and Z, or X, Y and Z periodicity

                    // X and Y
                    if (is_penultimate_periodic(0) && is_penultimate_periodic(1))
                    {
                        c_vector<int,DIM> ultimate_indices = grid_indices;
                        ultimate_indices(0) ++;
                        ultimate_indices(1) ++;

                        for (unsigned neighbour = 0 ; neighbour < neighbours.size() ; neighbour++)
                        {
                            if (IsBoxInDomain(ultimate_indices + neighbours[neighbour]))
                            {
                                local_boxes.insert(GetLinearIndex(ultimate_indices + neighbours[neighbour]));
                            }
                        }
                    }

                    // X and Z
                    if (is_penultimate_periodic(0) && is_penultimate_periodic(2))
                    {
                        c_vector<int,DIM> ultimate_indices = grid_indices;
                        ultimate_indices(0) ++;
                        ultimate_indices(2) ++;

                        for (unsigned neighbour = 0 ; neighbour < neighbours.size() ; neighbour++)
                        {
                            if (IsBoxInDomain(ultimate_indices + neighbours[neighbour]))
                            {
                                local_boxes.insert(GetLinearIndex(ultimate_indices + neighbours[neighbour]));
                            }
                        }
                    }

                    // Y and Z
                    if (is_penultimate_periodic(1) && is_penultimate_periodic(2))
                    {
                        c_vector<int,DIM> ultimate_indices = grid_indices;
                        ultimate_indices(1) ++;
                        ultimate_indices(2) ++;

                        for (unsigned neighbour = 0 ; neighbour < neighbours.size() ; neighbour++)
                        {
                            if (IsBoxInDomain(ultimate_indices + neighbours[neighbour]))
                            {
                                local_boxes.insert(GetLinearIndex(ultimate_indices + neighbours[neighbour]));
                            }
                        }
                    }

                    // X Y and Z
                    if (num_penultimate_periodic_dimensions == 3)
                    {
                        assert(mIsDomainPeriodic(0) && mIsDomainPeriodic(1) && mIsDomainPeriodic(1));
                        local_boxes.insert(0);
                    }

                    break;
                }
                default:
                {
                    NEVER_REACHED;
                }
            }
        }

        // Add the local boxes to the member vector
        mLocalBoxes.push_back(local_boxes);
    }
}

template<unsigned DIM>
std::set<unsigned> BoxCollection<DIM>::GetLocalBoxes(unsigned boxIndex)
{
    assert(boxIndex < mLocalBoxes.size());
    return mLocalBoxes[boxIndex];
}

template<unsigned DIM>
const c_vector<double, 2 * DIM>& BoxCollection<DIM>::rGetDomainSize() const
{
    return mDomainSize;
}

template<unsigned DIM>
void BoxCollection<DIM>::CalculateNodePairs(std::vector<Node<DIM>*>& rNodes,
                                            std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                            std::map<unsigned, std::set<unsigned> >& rNodeNeighbours)
{
    rNodePairs.clear();
    rNodeNeighbours.clear();

    // Ensure all boxes are empty
    EmptyBoxes();

    // Create an empty set of neighbours for each node, and add each node to its correct box
    for (unsigned node_index = 0; node_index < rNodes.size(); node_index++)
    {
        rNodeNeighbours[node_index] = std::set<unsigned>();

        unsigned box_index = CalculateContainingBox(rNodes[node_index]);
        mBoxes[box_index].AddNode(rNodes[node_index]);
    }

    for (unsigned i = 0; i < rNodes.size(); i++)
    {
        Node<DIM>* this_node = rNodes[i];
        unsigned node_index = this_node->GetIndex();

        // Get the box containing this node
        unsigned this_node_box_index = CalculateContainingBox(this_node);

        // Get the local boxes to this node
        std::set<unsigned> local_boxes_indices = GetLocalBoxes(this_node_box_index);

        // Loop over all the local boxes
        for (std::set<unsigned>::iterator box_iter = local_boxes_indices.begin(); box_iter != local_boxes_indices.end();
                box_iter++)
        {
            // Get the set of nodes contained in this box
            std::set<Node<DIM>*>& r_contained_nodes = mBoxes[*box_iter].rGetNodesContained();

            // Loop over these nodes
            for (typename std::set<Node<DIM>*>::iterator node_iter = r_contained_nodes.begin();
                    node_iter != r_contained_nodes.end(); ++node_iter)
            {
                // Get the index of the other node
                unsigned other_node_index = (*node_iter)->GetIndex();

                // If we're in the same box, then take care not to store the node pair twice
                if (*box_iter == this_node_box_index)
                {
                    if (other_node_index > this_node->GetIndex())
                    {
                        rNodePairs.push_back(std::pair<Node<DIM>*, Node<DIM>*>(this_node, (*node_iter)));
                        rNodeNeighbours[node_index].insert(other_node_index);
                        rNodeNeighbours[other_node_index].insert(node_index);
                    }
                }
                else
                {
                    rNodePairs.push_back(std::pair<Node<DIM>*, Node<DIM>*>(this_node, (*node_iter)));
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
