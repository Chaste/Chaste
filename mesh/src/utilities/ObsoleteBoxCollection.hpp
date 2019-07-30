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
#ifndef OBSOLETEBOXCOLLECTION_HPP_
#define OBSOLETEBOXCOLLECTION_HPP_

#include "Node.hpp"
#include "Element.hpp"
#include "Box.hpp"
#include <map>
#include <vector>


/**
 * A collection of 'boxes' partitioning the domain with information on which nodes are located in which box
 * Not archived - in cell_based constructed in NodeBasedCellPopulation constructor.
 *
 * This is to be merged with DistributedBoxCollection.
 */
template<unsigned DIM>
class ObsoleteBoxCollection
{
private:
    friend class TestObsoleteBoxCollection;

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

    /** Whether the domain is periodic */
    c_vector<bool, DIM> mIsDomainPeriodic;

    /** A fudge (box swelling) factor to deal with 32-bit floating point issues. */
    static const double msFudge;

    /**
     * @param gridIndices the grid indices (i), (i,j), or (i,j,k) depending on DIM
     * @return the linear index in row-major form
     */
    unsigned GetLinearIndex(c_vector<int, DIM> gridIndices);

    /**
     * @param linearIndex the linear index in row-major form
     * @return the grid indices (i), (i,j), or (i,j,k) depending on DIM
     */
    c_vector<unsigned, DIM> GetGridIndices(unsigned linearIndex);

    /**
     * @param gridIndices the coordinates (i), (i,j), or (i,j,k) depending on DIM
     * @return whether the box is in the domain
     */
    bool IsBoxInDomain(c_vector<unsigned, DIM> gridIndices);

    /**
     * This method is used for periodicity.  It is necessary to consider additional
     * boxes as being local to a candidate box if and only if the candidate box is
     * in the penultimate location of any dimension.
     *
     * @param gridIndices the coordinates (i), (i,j), or (i,j,k) depending on DIM
     * @return whether the box is in the penultimate location in each dimension
     */
    c_vector<bool,DIM> IsIndexPenultimate(c_vector<unsigned, DIM> gridIndices);

    /**
     * Helper function for SetupLocalBoxesHalfOnly() and SetupAllLocalBoxes().
     *
     * Accepts a vector of neighbours that should be considered for each box and populates
     * mLocalBoxes based on these neighbours.  The set of neighbours is either half or all
     * of the neighbours, depending on the function calling this method.
     *
     * @param rNeighbours a vector of neighbours
     */
    void SetupLocalBoxes(const std::vector<c_vector<int, DIM> >& rNeighbours);

public:

    /**
     * Constructor.
     *
     * @param boxWidth the width of each box (cut-off length in NodeBasedCellPopulation simulations)
     * @param domainSize the size of the domain, in the form (xmin, xmax, ymin, ymax) (etc)
     * @param isPeriodicInX whether the domain is periodic in the x direction
     * @param isPeriodicInY whether the domain is periodic in the y direction
     * @param isPeriodicInZ whether the domain is periodic in the z direction
     */
    ObsoleteBoxCollection(double boxWidth,
                          c_vector<double, 2*DIM> domainSize,
                          bool isPeriodicInX = false,
                          bool isPeriodicInY = false,
                          bool isPeriodicInZ = false);

    /**
     * Remove the list of nodes stored in each box.
     */
    void EmptyBoxes();

    /**
     * @return the index of the box this node is contained in.
     * @param pNode address of the node
     */
    unsigned CalculateContainingBox(Node<DIM>* pNode);

    /**
     * @return the index of the box a point is contained in
     * @param rLocation The point
     */
    unsigned CalculateContainingBox(c_vector<double, DIM>& rLocation);

    /**
     * @param boxIndex the index of the box to return
     * @return a reference to the box with global index boxIndex.
     */
    Box<DIM>& rGetBox(unsigned boxIndex);

    /**
     * @return the total (global) number of boxes.
     */
    unsigned GetNumBoxes();

    /**
     * Generates a list of vectors representing half the possible neighbour locations
     * in DIM-dimensions ({0} and {1} in 1D, rather than {-1}, {0}, {1}), and similar
     * in higher dimensions, and calls SetupLocalBoxes().
     */
    void SetupLocalBoxesHalfOnly();

    /**
     * Generates a list of vectors representing all the possible neighbour locations
     * in DIM-dimensions, and calls SetupLocalBoxes().
     */
    void SetupAllLocalBoxes();

    /**
     * Get the set of all the local boxes, i.e. itself and its nearest-neighbours.
     *
     * @param boxIndex the index of the box
     * @return the set containing the indices of boxes local to box boxIndex.  i.e. the box boxIndex itself and its nearest-neighbours.
     */
    std::set<unsigned>& rGetLocalBoxes(unsigned boxIndex);

    /**
     * Return the size of the domain partitioned by this collection
     *
     * @return the vector of the domain edges (xmin, xmax etc).
     */
    const c_vector<double, 2*DIM>& rGetDomainSize() const;

    /**
     * Compute all the pairs of (potentially) connected nodes for cell_based simulations, ie nodes which are in a
     * local box to the box containing the first node. **Note: the user still has to check that the node
     * pairs are less than the cut-off distance apart.** The pairs are checked so that index1 < index2,
     * so each connected pair of nodes is only in the set once.
     *
     * @param rNodes all the nodes to be considered
     * @param rNodePairs the return value, a set of pairs of nodes
     */
    void CalculateNodePairs(std::vector<Node<DIM>*>& rNodes,
                            std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);

    /**
     *  Helper function for transition to distributed version of box collection.
     *
     *  This method allows us to insert of dummy test for whether a box/point/node/location is located
     *  on the local process.
     *
     * @param pNode the node to test.
     * @return true (whether the point at pNode->rGetLocation() is owned on this process)
     */
    bool IsOwned(Node<DIM>* pNode)
    {
        return true;
    }

    /**
     *  Helper function for transition to distributed version of box collection.
     *
     *  This method allows us to insert of dummy test for whether a box/point/node/location is located
     *  on the local process.
     * @param globalIndex the global index of the box.
     * @return true (whether the point at pNode->rGetLocation() is owned on this process)
     */
    bool IsBoxOwned(unsigned globalIndex)
    {
        return true;
    }
};

#endif /*OBSOLETEBOXCOLLECTION_HPP_*/
