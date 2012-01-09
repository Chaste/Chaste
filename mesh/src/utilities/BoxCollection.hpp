/*

Copyright (C) University of Oxford, 2005-2012

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
#ifndef BOXCOLLECTION_HPP_
#define BOXCOLLECTION_HPP_

#include "Node.hpp"
#include "Element.hpp"


/**
 * A small class for a nD 'box' defined by its min/max x/y/z values which
 * can contains a list of nodes and elements located in that box
 */
template<unsigned DIM>
class Box
{
private:
    /** Coordinates of the box, in the form (for 2D) (xmin, xmax, ymin, ymax) (etc). */
    c_vector<double, 2*DIM> mMinAndMaxValues;

    /** Nodes contained in this box. */
    std::set< Node<DIM>* > mNodesContained;

    /** Elements contained in this box. */
    std::set< Element<DIM,DIM>* > mElementsContained;

public:

    /**
     * Constructor just takes in the extremal values of the box.
     *
     * @param rMinAndMaxValues the extremal values. Of the from (for 2D, etc): xmin, xmax, ymin, ymax
     */
    Box(c_vector<double, 2*DIM>& rMinAndMaxValues);

    /** Get the coordinates of the box, in the form (for 2D) (xmin, xmax, ymin, ymax) (etc). */
    c_vector<double, 2*DIM>& rGetMinAndMaxValues();

    /**
     * Add a node to this box.
     * @param pNode address of the node to be added
     */
    void AddNode(Node<DIM>* pNode);

    /**
     * Remove a node from this box.
     * @param pNode address of the node to be removed
     */
    void RemoveNode(Node<DIM>* pNode);

    /**
     * An element to this box.
     * @param pElement address of the element to be added
     */
    void AddElement(Element<DIM,DIM>* pElement);

    /** Get all the nodes in this box. */
    std::set< Node<DIM>* >& rGetNodesContained();

    /** Get all the elements in this box. */
    std::set< Element<DIM,DIM>* >& rGetElementsContained();
};


/**
 * A collection of 'boxes' partitioning the domain with information on which nodes are located in which box
 * Not archived - in cell_based constructed in NodeBasedCellPopulation constructor.
 */
template<unsigned DIM>
class BoxCollection
{
private:
    /** A vector of boxes to store rough node/element positions. */
    std::vector< Box<DIM> > mBoxes;

    /** The domain being partitioned. */
    c_vector<double, 2*DIM> mDomainSize;

    /** The width of each box. */
    double mBoxWidth;

    /** Number of boxes in each direction. */
    c_vector<unsigned, DIM> mNumBoxesEachDirection;

    /** The boxes local (itself and nearest neighbour) to a given box. */
    std::vector< std::set<unsigned> > mLocalBoxes;

public:

    /**
     * Constructor
     *
     * @param boxWidth the width of each box (cut-off length in NodeBasedCellPopulation simulations)
     * @param domainSize the size of the domain, in the form (xmin, xmax, ymin, ymax) (etc)
     */
    BoxCollection(double boxWidth, c_vector<double, 2*DIM> domainSize);

    /**
     * Calculate which box this node is contained in.
     * @param pNode address of the node
     */
    unsigned CalculateContainingBox(Node<DIM>* pNode);

    /**
     * Calculate which box a point is contained in
     * @param rLocation The point
     */
    unsigned CalculateContainingBox(c_vector<double, DIM>& rLocation);

    /**
     * Get a box.
     * @param boxIndex the index of the box to return
     */
    Box<DIM>& rGetBox(unsigned boxIndex);

    /** Get the number of boxes. */
    unsigned GetNumBoxes();

    /** Set up the local boxes (ie itself and its nearest-neighbours) for each of the boxes.
     *  Just set up half of the local boxes (for example, in 1D, local boxes for box0 = {1}
     *  local boxes for box1 = {2} not {0,2}, and so on. Similar to 2d, 3d.
     */
    void SetupLocalBoxesHalfOnly();

    /** Set up the local boxes (ie itself and its nearest-neighbours) for each of the boxes. */
    void SetupAllLocalBoxes();

    /**
     * Get the set of all the local boxes, i.e. itself and its nearest-neighbours.
     * @param boxIndex the index of the box
     */
    std::set<unsigned> GetLocalBoxes(unsigned boxIndex);

    /**
     *  Compute all the pairs of (potentially) connected nodes for cell_based simulations, ie nodes which are in a
     *  local box to the box containing the first node. **Note: the user still has to check that the node
     *  pairs are less than the cut-off distance apart.** The pairs are checked so that index1 < index2,
     *  so each connected pair of nodes is only in the set once.
     *
     *  @param rNodes all the nodes to be consider
     *  @param rNodePairs the return value, a set of pairs of nodes
     */
    void CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::set<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);
};


#endif /*BOXCOLLECTION_HPP_*/
