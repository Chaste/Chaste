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

#include "PottsMeshGenerator.hpp"

#include <boost/scoped_array.hpp>

template<unsigned DIM>
PottsMeshGenerator<DIM>::PottsMeshGenerator(unsigned numNodesAcross, unsigned numElementsAcross, unsigned elementWidth,
                                            unsigned numNodesUp, unsigned numElementsUp, unsigned elementHeight,
                                            unsigned numNodesDeep, unsigned numElementsDeep, unsigned elementDepth,
                                            bool startAtBottomLeft, bool isPeriodicInX, bool isPeriodicInY ,bool isPeriodicInZ)
{
    assert(numNodesAcross > 0);
    assert(numNodesUp > 0);
    assert(numNodesDeep > 0);

    assert(numElementsAcross*elementWidth <= numNodesAcross);
    assert(numElementsUp*elementHeight <= numNodesUp);
    assert(numElementsDeep*elementDepth <= numNodesDeep);

    std::vector<Node<DIM>*> nodes;
    std::vector<PottsElement<DIM>*>  elements;
    std::vector<std::set<unsigned> > moore_neighbours;
    std::vector<std::set<unsigned> > von_neumann_neighbours;

    unsigned num_nodes = numNodesAcross*numNodesUp*numNodesDeep;

    unsigned next_node_index = 0;
    boost::scoped_array<unsigned> node_indices(new unsigned[elementWidth*elementHeight*elementDepth]);
    unsigned element_index;

    unsigned index_offset = 0;

    if (!startAtBottomLeft) // Elements in centre of mesh
    {
        // Calculate the width of the medium on the edge and offset the node index so that the elements are in the centre of the mesh.
        unsigned across_gap = (numNodesAcross -  numElementsAcross*elementWidth)/2;
        unsigned up_gap = (numNodesUp -  numElementsUp*elementHeight)/2;
        unsigned deep_gap = (numNodesDeep -  numElementsDeep*elementDepth)/2;

        index_offset = deep_gap*numNodesAcross*numNodesUp + up_gap*numNodesAcross + across_gap;
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
                    is_boundary_node = (j==0 || j==numNodesUp-1 || (i==0 && !isPeriodicInX) || (i==numNodesAcross-1 && !isPeriodicInX) ) ? true : false;
                }
                if (DIM==3)
                {
                    is_boundary_node = (j==0 || j==numNodesUp-1 || (i==0 && !isPeriodicInX) || (i==numNodesAcross-1 && !isPeriodicInX) || k==0 || k==numNodesDeep-1) ? true : false;
                }
                Node<DIM>* p_node = new Node<DIM>(next_node_index, is_boundary_node, i, j, k);
                nodes.push_back(p_node);
                next_node_index++;
            }
        }
    }
    assert(nodes.size()==num_nodes);

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

    /*
     * Create the neighborhoods of each node
     */

    moore_neighbours.resize(num_nodes);
    von_neumann_neighbours.resize(num_nodes);

    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        // Clear the set of neighbouring node indices

        moore_neighbours[node_index].clear();

        switch (DIM)
        {
            case 2:
            {
                assert(DIM == 2);
                /*
                 * This stores the available neighbours using the following numbering:
                 *
                 *  1-----0-----7
                 *  |     |     |
                 *  |     |     |
                 *  2-----x-----6
                 *  |     |     |
                 *  |     |     |
                 *  3-----4-----5
                 */

                // Create a vector of possible neighbouring node indices
                std::vector<unsigned> moore_neighbour_indices_vector(8, node_index);
                moore_neighbour_indices_vector[0] += numNodesAcross;
                moore_neighbour_indices_vector[1] += numNodesAcross - 1;
                moore_neighbour_indices_vector[2] -= 1;
                moore_neighbour_indices_vector[3] -= numNodesAcross + 1;
                moore_neighbour_indices_vector[4] -= numNodesAcross;
                moore_neighbour_indices_vector[5] -= numNodesAcross - 1;
                moore_neighbour_indices_vector[6] += 1;
                moore_neighbour_indices_vector[7] += numNodesAcross + 1;

                // Work out whether this node lies on any edge of the mesh
                bool on_south_edge = (node_index < numNodesAcross);
                bool on_north_edge = ((int) node_index > (int)numNodesAcross*((int)numNodesUp - 1) - 1);
                bool on_west_edge = (node_index%numNodesAcross == 0);
                bool on_east_edge = (node_index%numNodesAcross == numNodesAcross - 1);

                if (isPeriodicInX)
                {
                    if (on_west_edge)
                    {
                        moore_neighbour_indices_vector[1] = node_index + 2*numNodesAcross - 1;
                        moore_neighbour_indices_vector[2] = node_index + numNodesAcross - 1;
                        moore_neighbour_indices_vector[3] = node_index - 1;
                    }

                    if (on_east_edge)
                    {
                        moore_neighbour_indices_vector[5] = node_index - 2*numNodesAcross + 1;
                        moore_neighbour_indices_vector[6] = node_index - numNodesAcross + 1;
                        moore_neighbour_indices_vector[7] = node_index + 1;
                    }
                }

                if (isPeriodicInY)
                {
                    if (on_north_edge)
                    {
                        moore_neighbour_indices_vector[0] = node_index - numNodesAcross*(numNodesUp-1);
                        moore_neighbour_indices_vector[1] = moore_neighbour_indices_vector[0] - 1;
                        moore_neighbour_indices_vector[7] = moore_neighbour_indices_vector[0] + 1;

                        if (on_west_edge)
                        {
                            moore_neighbour_indices_vector[1] = moore_neighbour_indices_vector[1] + numNodesAcross;
                        }
                        if (on_east_edge)
                        {
                            moore_neighbour_indices_vector[7] = moore_neighbour_indices_vector[7] - numNodesAcross;
                        }
                    }

                    if (on_south_edge)
                    {
                        moore_neighbour_indices_vector[4] = node_index + numNodesAcross*(numNodesUp-1);
                        moore_neighbour_indices_vector[3] = moore_neighbour_indices_vector[4] - 1;
                        moore_neighbour_indices_vector[5] = moore_neighbour_indices_vector[4] + 1;

                        if (on_west_edge)
                        {
                            moore_neighbour_indices_vector[3] = moore_neighbour_indices_vector[3] + numNodesAcross;
                        }
                        if (on_east_edge)
                        {
                            moore_neighbour_indices_vector[5] = moore_neighbour_indices_vector[5] - numNodesAcross;
                        }
                    }
                }

                if (isPeriodicInX)
                {
                    on_east_edge = false;
                    on_west_edge = false;
                }
                if (isPeriodicInY)
                {
                    on_south_edge = false;
                    on_north_edge = false;
                }

                // Create a vector of booleans for which neighbours are available
                // Use the order N, NW, W, SW, S, SE, E, NE
                std::vector<bool> available_neighbours = std::vector<bool>(8, true);
                available_neighbours[0] = !on_north_edge;
                available_neighbours[1] = !(on_north_edge || on_west_edge);
                available_neighbours[2] = !on_west_edge;
                available_neighbours[3] = !(on_south_edge || on_west_edge);
                available_neighbours[4] = !on_south_edge;
                available_neighbours[5] = !(on_south_edge || on_east_edge);
                available_neighbours[6] = !on_east_edge;
                available_neighbours[7] = !(on_north_edge || on_east_edge);

                // Using neighbour_indices_vector and available_neighbours, store the indices of all available neighbours to the set all_neighbours
                for (unsigned i=0; i<8; i++)
                {
                    if (available_neighbours[i])
                    {
                        assert(moore_neighbour_indices_vector[i] < nodes.size());
                        moore_neighbours[node_index].insert(moore_neighbour_indices_vector[i]);
                    }
                }
                break;
            }
            case 3:
            {
                assert(DIM ==3);
                /*
                 * This stores the available neighbours using the following numbering:
                 *                      FRONT           BACK
                 *  1-----0-----7   10-----9---16   19----18---25
                 *  |     |     |   |     |     |   |     |     |
                 *  |     |     |   |     |     |   |     |     |
                 *  2-----x-----6   11----8-----15  20----17---24
                 *  |     |     |   |     |     |   |     |     |
                 *  |     |     |   |     |     |   |     |     |
                 *  3-----4-----5   12----13----14  21---22----23
                 */

                // Create a vector of possible neighbouring node indices
                std::vector<unsigned> moore_neighbour_indices_vector(26, node_index);
                moore_neighbour_indices_vector[0] += numNodesAcross;
                moore_neighbour_indices_vector[1] += numNodesAcross - 1;
                moore_neighbour_indices_vector[2] -= 1;
                moore_neighbour_indices_vector[3] -= numNodesAcross + 1;
                moore_neighbour_indices_vector[4] -= numNodesAcross;
                moore_neighbour_indices_vector[5] -= numNodesAcross - 1;
                moore_neighbour_indices_vector[6] += 1;
                moore_neighbour_indices_vector[7] += numNodesAcross + 1;
                moore_neighbour_indices_vector[8] -= numNodesAcross*numNodesUp;
                for (unsigned i=9; i<17; i++)
                {
                    moore_neighbour_indices_vector[i] = moore_neighbour_indices_vector[i-9]-numNodesAcross*numNodesUp;
                }
                moore_neighbour_indices_vector[17] += numNodesAcross*numNodesUp;
                for (unsigned i=18; i<26; i++)
                {
                    moore_neighbour_indices_vector[i]=moore_neighbour_indices_vector[i-18]+numNodesAcross*numNodesUp;
                }

                // Work out whether this node lies on any edge of the mesh
                bool on_south_edge = (node_index%(numNodesAcross*numNodesUp)<numNodesAcross);
                bool on_north_edge = (node_index%(numNodesAcross*numNodesUp)>(numNodesAcross*numNodesUp-numNodesAcross-1));
                bool on_west_edge = (node_index%numNodesAcross == 0);
                bool on_east_edge = (node_index%numNodesAcross == numNodesAcross - 1);
                bool on_front_edge = (node_index < numNodesAcross*numNodesUp);
                bool on_back_edge = (node_index > numNodesAcross*numNodesUp*(numNodesDeep-1)-1);

                if (isPeriodicInX)
                {
                    if (on_west_edge)
                    {
                        moore_neighbour_indices_vector[1] = node_index + 2*numNodesAcross - 1;
                        moore_neighbour_indices_vector[2] = node_index + numNodesAcross - 1;
                        moore_neighbour_indices_vector[3] = node_index - 1;

                        moore_neighbour_indices_vector[10] = moore_neighbour_indices_vector[1] - numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[11] = moore_neighbour_indices_vector[2] - numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[12] = moore_neighbour_indices_vector[3] - numNodesAcross*numNodesUp;

                        moore_neighbour_indices_vector[19] = moore_neighbour_indices_vector[1] + numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[20] = moore_neighbour_indices_vector[2] + numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[21] = moore_neighbour_indices_vector[3] + numNodesAcross*numNodesUp;
                    }

                    if (on_east_edge)
                    {
                        moore_neighbour_indices_vector[5] = node_index - 2*numNodesAcross + 1;
                        moore_neighbour_indices_vector[6] = node_index - numNodesAcross + 1;
                        moore_neighbour_indices_vector[7] = node_index + 1;

                        moore_neighbour_indices_vector[14] = moore_neighbour_indices_vector[5] - numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[15] = moore_neighbour_indices_vector[6] - numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[16] = moore_neighbour_indices_vector[7] - numNodesAcross*numNodesUp;

                        moore_neighbour_indices_vector[23] = moore_neighbour_indices_vector[5] + numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[24] = moore_neighbour_indices_vector[6] + numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[25] = moore_neighbour_indices_vector[7] + numNodesAcross*numNodesUp;
                    }
                }

                if (isPeriodicInY)
                {
                    if (on_north_edge)
                    {
                        moore_neighbour_indices_vector[0] = node_index - numNodesAcross*(numNodesUp-1);
                        moore_neighbour_indices_vector[1] = moore_neighbour_indices_vector[0] - 1;
                        moore_neighbour_indices_vector[7] = moore_neighbour_indices_vector[0] + 1;

                        moore_neighbour_indices_vector[10] = moore_neighbour_indices_vector[1] - numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[9] = moore_neighbour_indices_vector[0] - numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[16] = moore_neighbour_indices_vector[7] - numNodesAcross*numNodesUp;

                        moore_neighbour_indices_vector[19] = moore_neighbour_indices_vector[1] + numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[18] = moore_neighbour_indices_vector[0] + numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[25] = moore_neighbour_indices_vector[7] + numNodesAcross*numNodesUp;

                        if (on_west_edge)
                        {
                            moore_neighbour_indices_vector[1] = moore_neighbour_indices_vector[1] + numNodesAcross;
                            moore_neighbour_indices_vector[10] = moore_neighbour_indices_vector[10] + numNodesAcross;
                            moore_neighbour_indices_vector[19] = moore_neighbour_indices_vector[19] + numNodesAcross;
                        }
                        if (on_east_edge)
                        {
                            moore_neighbour_indices_vector[7] = moore_neighbour_indices_vector[7] - numNodesAcross;
                            moore_neighbour_indices_vector[16] = moore_neighbour_indices_vector[16] - numNodesAcross;
                            moore_neighbour_indices_vector[25] = moore_neighbour_indices_vector[25] - numNodesAcross;
                        }
                    }

                    if (on_south_edge)
                    {
                        moore_neighbour_indices_vector[4] = node_index + numNodesAcross*(numNodesUp-1);
                        moore_neighbour_indices_vector[3] = moore_neighbour_indices_vector[4] - 1;
                        moore_neighbour_indices_vector[5] = moore_neighbour_indices_vector[4] + 1;

                        moore_neighbour_indices_vector[12] = moore_neighbour_indices_vector[3] - numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[13] = moore_neighbour_indices_vector[4] - numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[14] = moore_neighbour_indices_vector[5] - numNodesAcross*numNodesUp;

                        moore_neighbour_indices_vector[21] = moore_neighbour_indices_vector[3] + numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[22] = moore_neighbour_indices_vector[4] + numNodesAcross*numNodesUp;
                        moore_neighbour_indices_vector[23] = moore_neighbour_indices_vector[5] + numNodesAcross*numNodesUp;

                        if (on_west_edge)
                        {
                            moore_neighbour_indices_vector[3] = moore_neighbour_indices_vector[3] + numNodesAcross;
                            moore_neighbour_indices_vector[12] = moore_neighbour_indices_vector[12] + numNodesAcross;
                            moore_neighbour_indices_vector[21] = moore_neighbour_indices_vector[21] + numNodesAcross;
                        }
                        if (on_east_edge)
                        {
                            moore_neighbour_indices_vector[5] = moore_neighbour_indices_vector[5] - numNodesAcross;
                            moore_neighbour_indices_vector[14] = moore_neighbour_indices_vector[14] - numNodesAcross;
                            moore_neighbour_indices_vector[23] = moore_neighbour_indices_vector[23] - numNodesAcross;
                        }
                    }
                }

                if (isPeriodicInZ)
                {
                    if (on_back_edge)
                    {
                        moore_neighbour_indices_vector[17] = node_index - numNodesAcross*numNodesUp*(numNodesDeep-1);
                        moore_neighbour_indices_vector[20] = moore_neighbour_indices_vector[17] - 1;
                        moore_neighbour_indices_vector[24] = moore_neighbour_indices_vector[17] + 1;

                        moore_neighbour_indices_vector[21] = moore_neighbour_indices_vector[20] - numNodesAcross;
                        moore_neighbour_indices_vector[22] = moore_neighbour_indices_vector[17] - numNodesAcross;
                        moore_neighbour_indices_vector[23] = moore_neighbour_indices_vector[24] - numNodesAcross;

                        moore_neighbour_indices_vector[19] = moore_neighbour_indices_vector[20] + numNodesAcross;
                        moore_neighbour_indices_vector[18] = moore_neighbour_indices_vector[17] + numNodesAcross;
                        moore_neighbour_indices_vector[25] = moore_neighbour_indices_vector[24] + numNodesAcross;

                        if (on_west_edge)
                        {
                            moore_neighbour_indices_vector[19] = moore_neighbour_indices_vector[19] + numNodesAcross;
                            moore_neighbour_indices_vector[20] = moore_neighbour_indices_vector[20] + numNodesAcross;
                            moore_neighbour_indices_vector[21] = moore_neighbour_indices_vector[21] + numNodesAcross;
                        }
                        if (on_east_edge)
                        {
                            moore_neighbour_indices_vector[23] = moore_neighbour_indices_vector[23] - numNodesAcross;
                            moore_neighbour_indices_vector[24] = moore_neighbour_indices_vector[24] - numNodesAcross;
                            moore_neighbour_indices_vector[25] = moore_neighbour_indices_vector[25] - numNodesAcross;
                        }

                        if (on_south_edge)
                        {
                            moore_neighbour_indices_vector[21] = moore_neighbour_indices_vector[21] + numNodesAcross*numNodesUp;
                            moore_neighbour_indices_vector[22] = moore_neighbour_indices_vector[22] + numNodesAcross*numNodesUp;
                            moore_neighbour_indices_vector[23] = moore_neighbour_indices_vector[23] + numNodesAcross*numNodesUp;
                        }

                        if (on_north_edge)
                        {
                            moore_neighbour_indices_vector[18] = moore_neighbour_indices_vector[18] - numNodesAcross*numNodesUp;
                            moore_neighbour_indices_vector[19] = moore_neighbour_indices_vector[19] - numNodesAcross*numNodesUp;
                            moore_neighbour_indices_vector[25] = moore_neighbour_indices_vector[25] - numNodesAcross*numNodesUp;
                        }
                    }

                    if (on_front_edge)
                    {
                        moore_neighbour_indices_vector[8] = node_index + numNodesAcross*numNodesUp*(numNodesDeep-1);
                        moore_neighbour_indices_vector[11] = moore_neighbour_indices_vector[8] - 1;
                        moore_neighbour_indices_vector[15] = moore_neighbour_indices_vector[8] + 1;

                        moore_neighbour_indices_vector[12] = moore_neighbour_indices_vector[11] - numNodesAcross;
                        moore_neighbour_indices_vector[13] = moore_neighbour_indices_vector[8] - numNodesAcross;
                        moore_neighbour_indices_vector[14] = moore_neighbour_indices_vector[15] - numNodesAcross;

                        moore_neighbour_indices_vector[10] = moore_neighbour_indices_vector[11] + numNodesAcross;
                        moore_neighbour_indices_vector[9] = moore_neighbour_indices_vector[8] + numNodesAcross;
                        moore_neighbour_indices_vector[16] = moore_neighbour_indices_vector[15] + numNodesAcross;

                        if (on_west_edge)
                        {
                            moore_neighbour_indices_vector[10] = moore_neighbour_indices_vector[10] + numNodesAcross;
                            moore_neighbour_indices_vector[11] = moore_neighbour_indices_vector[11] + numNodesAcross;
                            moore_neighbour_indices_vector[12] = moore_neighbour_indices_vector[12] + numNodesAcross;
                        }
                        if (on_east_edge)
                        {
                            moore_neighbour_indices_vector[14] = moore_neighbour_indices_vector[14] - numNodesAcross;
                            moore_neighbour_indices_vector[15] = moore_neighbour_indices_vector[15] - numNodesAcross;
                            moore_neighbour_indices_vector[16] = moore_neighbour_indices_vector[16] - numNodesAcross;
                        }

                        if (on_south_edge)
                        {
                            moore_neighbour_indices_vector[12] = moore_neighbour_indices_vector[12] + numNodesAcross*numNodesUp;
                            moore_neighbour_indices_vector[13] = moore_neighbour_indices_vector[13] + numNodesAcross*numNodesUp;
                            moore_neighbour_indices_vector[14] = moore_neighbour_indices_vector[14] + numNodesAcross*numNodesUp;
                        }

                        if (on_north_edge)
                        {
                            moore_neighbour_indices_vector[9] = moore_neighbour_indices_vector[9] - numNodesAcross*numNodesUp;
                            moore_neighbour_indices_vector[10] = moore_neighbour_indices_vector[10] - numNodesAcross*numNodesUp;
                            moore_neighbour_indices_vector[16] = moore_neighbour_indices_vector[16] - numNodesAcross*numNodesUp;
                        }
                    }
                }

                if (isPeriodicInX)
                {
                    on_east_edge = false;
                    on_west_edge = false;
                }
                if (isPeriodicInY)
                {
                    on_south_edge = false;
                    on_north_edge = false;
                }
                if (isPeriodicInZ)
                {
                    on_front_edge = false;
                    on_back_edge = false;
                }

                // Create a vector of booleans for which neighbours are available
                // Use the order N, NW, W, SW, S, SE, E, NE
                std::vector<bool> available_neighbours = std::vector<bool>(26, true);
                available_neighbours[0] = !on_north_edge;
                available_neighbours[1] = !(on_north_edge || on_west_edge);
                available_neighbours[2] = !on_west_edge;
                available_neighbours[3] = !(on_south_edge || on_west_edge);
                available_neighbours[4] = !on_south_edge;
                available_neighbours[5] = !(on_south_edge || on_east_edge);
                available_neighbours[6] = !on_east_edge;
                available_neighbours[7] = !(on_north_edge || on_east_edge);
                available_neighbours[8] = !(on_front_edge);
                for (unsigned i=9; i<17; i++)
                {
                    available_neighbours[i] = (available_neighbours[i-9] && !(on_front_edge));
                }
                available_neighbours[17] = !(on_back_edge);
                for (unsigned i=18; i<26; i++)
                {
                    available_neighbours[i] = (available_neighbours[i-18] && !(on_back_edge));
                }

                // Using neighbour_indices_vector and available_neighbours, store the indices of all available neighbours to the set all_neighbours
                for (unsigned i=0; i<26; i++)
                {
                    if (available_neighbours[i] && moore_neighbour_indices_vector[i] < numNodesAcross*numNodesUp*numNodesDeep)
                    {
                        assert(moore_neighbour_indices_vector[i] < nodes.size());
                        moore_neighbours[node_index].insert(moore_neighbour_indices_vector[i]);
                    }
                }
                break;
            }
            default:
                NEVER_REACHED;
        }

        // Clear the set of neighbouring node indices
        von_neumann_neighbours[node_index].clear();

        switch (DIM)
        {
            case 2:
            {
                assert(DIM == 2);
                /*
                 * This stores the available neighbours using the following numbering:
                 *
                 *        0
                 *        |
                 *        |
                 *  1-----x-----3
                 *        |
                 *        |
                 *        2
                 */
                // Create a vector of possible neighbouring node indices
                std::vector<unsigned> von_neumann_neighbour_indices_vector(4, node_index);
                von_neumann_neighbour_indices_vector[0] += numNodesAcross;
                von_neumann_neighbour_indices_vector[1] -= 1;
                von_neumann_neighbour_indices_vector[2] -= numNodesAcross;
                von_neumann_neighbour_indices_vector[3] += 1;

                // Work out whether this node lies on any edge of the mesh
                bool on_south_edge = (node_index < numNodesAcross);
                bool on_north_edge = ((int)node_index > (int)numNodesAcross*((int)numNodesUp - 1) - 1);
                bool on_west_edge = (node_index%numNodesAcross == 0);
                bool on_east_edge = (node_index%numNodesAcross == numNodesAcross - 1);

                if (isPeriodicInX)
                {
                    if (on_west_edge)
                    {
                        von_neumann_neighbour_indices_vector[1] = node_index + numNodesAcross - 1;
                        on_west_edge = false;
                    }

                    if (on_east_edge)
                    {
                        von_neumann_neighbour_indices_vector[3] = node_index - numNodesAcross + 1;
                        on_east_edge = false;
                    }
                }

                if (isPeriodicInY)
                {
                    if (on_north_edge)
                    {
                        von_neumann_neighbour_indices_vector[0] = node_index - numNodesAcross*(numNodesUp-1);
                        on_north_edge = false;
                    }

                    if (on_south_edge)
                    {
                        von_neumann_neighbour_indices_vector[2] = node_index + numNodesAcross*(numNodesUp-1);
                        on_south_edge = false;
                    }
                }

                // Create a vector of booleans for which neighbours are available
                // Use the order N, W, S, E
                std::vector<bool> available_neighbours = std::vector<bool>(2*DIM, true);
                available_neighbours[0] = !on_north_edge;
                available_neighbours[1] = !on_west_edge;
                available_neighbours[2] = !on_south_edge;
                available_neighbours[3] = !on_east_edge;

                // Using von_neumann_neighbour_indices_vector and available_neighbours, store the indices of all available neighbours to the set all_neighbours
                for (unsigned i=0; i<4; i++)
                {
                    if (available_neighbours[i])
                    {
                        assert(von_neumann_neighbour_indices_vector[i] < nodes.size());
                        von_neumann_neighbours[node_index].insert(von_neumann_neighbour_indices_vector[i]);
                    }
                }
                break;
            }
            case 3:
            {
                assert(DIM == 3);

                /*
                 * This stores the available neighbours using the following numbering:
                 *
                 *        0  5
                 *        | /
                 *        |/
                 *  1-----x-----3
                 *      / |
                 *     /  |
                 *    4   2
                 */
                // Create a vector of possible neighbouring node indices
                std::vector<unsigned> von_neumann_neighbour_indices_vector(6, node_index);
                von_neumann_neighbour_indices_vector[0] += numNodesAcross;
                von_neumann_neighbour_indices_vector[1] -= 1;
                von_neumann_neighbour_indices_vector[2] -= numNodesAcross;
                von_neumann_neighbour_indices_vector[3] += 1;
                von_neumann_neighbour_indices_vector[4] -= numNodesAcross*numNodesUp;
                von_neumann_neighbour_indices_vector[5] += numNodesAcross*numNodesUp;

                // Work out whether this node lies on any edge of the mesh
                bool on_south_edge = (node_index%(numNodesAcross*numNodesUp)<numNodesAcross);
                bool on_north_edge = (node_index%(numNodesAcross*numNodesUp)>(numNodesAcross*numNodesUp-numNodesAcross-1));
                bool on_west_edge = (node_index%numNodesAcross== 0);
                bool on_east_edge = (node_index%numNodesAcross == numNodesAcross - 1);
                bool on_front_edge = (node_index < numNodesAcross*numNodesUp);
                bool on_back_edge = (node_index > numNodesAcross*numNodesUp*(numNodesDeep-1)-1);

                if (isPeriodicInX)
                {
                    if (on_west_edge)
                    {
                        von_neumann_neighbour_indices_vector[1] = node_index + numNodesAcross - 1;
                        on_west_edge = false;
                    }

                    if (on_east_edge)
                    {
                        von_neumann_neighbour_indices_vector[3] = node_index - numNodesAcross + 1;
                        on_east_edge = false;
                    }
                }

                if (isPeriodicInY)
                {
                    if (on_north_edge)
                    {
                        von_neumann_neighbour_indices_vector[0] = node_index - numNodesAcross*(numNodesUp-1);
                        on_north_edge = false;
                    }

                    if (on_south_edge)
                    {
                        von_neumann_neighbour_indices_vector[2] = node_index + numNodesAcross*(numNodesUp-1);
                        on_south_edge = false;
                    }
                }

                if (isPeriodicInZ)
                {
                    if (on_front_edge)
                    {
                        von_neumann_neighbour_indices_vector[4] = node_index + numNodesAcross*numNodesUp*(numNodesDeep-1);
                        on_front_edge = false;
                    }

                    if (on_back_edge)
                    {
                        von_neumann_neighbour_indices_vector[5] = node_index - numNodesAcross*numNodesUp*(numNodesDeep-1);
                        on_back_edge = false;
                    }
                }

                // Create a vector of booleans for which neighbours are available
                // Use the order N, W, S, E, F, B
                std::vector<bool> available_neighbours = std::vector<bool>(2*DIM, true);
                available_neighbours[0] = !on_north_edge;
                available_neighbours[1] = !on_west_edge;
                available_neighbours[2] = !on_south_edge;
                available_neighbours[3] = !on_east_edge;
                available_neighbours[4] = !on_front_edge;
                available_neighbours[5] = !on_back_edge;

                // Using von_neumann_neighbour_indices_vector and available_neighbours, store the indices of all available neighbours to the set all_neighbours
                for (unsigned i=0; i<6; i++)
                {
                    if (available_neighbours[i] && von_neumann_neighbour_indices_vector[i]<numNodesAcross*numNodesUp*numNodesDeep)
                    {
                        assert(von_neumann_neighbour_indices_vector[i] < nodes.size());
                        von_neumann_neighbours[node_index].insert(von_neumann_neighbour_indices_vector[i]);
                    }
                }
                break;
            }
            default:
                NEVER_REACHED;
        }
    }

    mpMesh = new PottsMesh<DIM>(nodes, elements, von_neumann_neighbours, moore_neighbours);
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

// Explicit instantiation
template class PottsMeshGenerator<1>;
template class PottsMeshGenerator<2>;
template class PottsMeshGenerator<3>;
