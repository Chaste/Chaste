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
#ifndef DISTRIBUTEDBOXCOLLECTION_HPP_
#define DISTRIBUTEDBOXCOLLECTION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>

#include "Node.hpp"
#include "Element.hpp"
#include "Box.hpp"
#include "PetscTools.hpp"
#include "DistributedVectorFactory.hpp"
#include <map>
#include <vector>


/**
 * A collection of 'boxes' partitioning the domain with information on which nodes are located in which box.
 */
template<unsigned DIM>
class DistributedBoxCollection
{
private:
    friend class TestDistributedBoxCollection;

    /** A vector of boxes to store rough node/element positions. */
    std::vector< Box<DIM> > mBoxes;

    /** A vector of boxes owned on other processes sharing a boundary with this process */
    std::vector< Box<DIM> > mHaloBoxes;

    /** A vector of the global indices of boxes, owned by this process, but on a boundary with right process. */
    std::vector<unsigned> mHalosRight;

    /** A vector of the global indices of boxes, owned by this process, but on a boundary with left process. */
    std::vector<unsigned> mHalosLeft;

    /** Set of Nodes that are halos of adjacent right process, but lie locally */
    std::vector<unsigned> mHaloNodesRight;

    /** Set of Nodes that are halos of adjacent left process, but lie locally */
    std::vector<unsigned> mHaloNodesLeft;

    /** Map of global to local indices of halo boxes in mHaloBoxes. **/
    std::map<unsigned, unsigned> mHaloBoxesMapping;

    /** The domain being partitioned. */
    c_vector<double, 2*DIM> mDomainSize;

    /** The width of each box. */
    double mBoxWidth;

    /** The total number of boxes across all processes */
    unsigned mNumBoxes;

    /** Number of boxes in each direction. */
    c_vector<unsigned, DIM> mNumBoxesEachDirection;

    /** Number of boxes in a face (1 in 1d, mNumBoxesEachDirection(0) in 2d, mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1) in 3d */
    unsigned mNumBoxesInAFace;

    /** The boxes local (itself and nearest neighbour) to a given box. */
    std::vector< std::set<unsigned> > mLocalBoxes;

    /** The smallest index of the boxes owned by this process. */
    unsigned mMinBoxIndex;

    /** The largest index of the boxes owned by this process. */
    unsigned mMaxBoxIndex;

    /** Whether the domain is periodic in the X dimension Note this currently only works for DIM=2.*/
    bool mIsPeriodicInX;

    /** Whether the local boxes have been setup or not. */
    bool mAreLocalBoxesSet;

    /** A fudge (box swelling) factor to deal with 32-bit floating point issues. */
    static const double msFudge;

    /** A distributed vector factory that governs ownership of rows of boxes */
    DistributedVectorFactory* mpDistributedBoxStackFactory;

    /** A flag that can be set to not save rNodeNeighbours in CalculateNodePairs - for efficiency */
    bool mCalculateNodeNeighbours;

    /**
     * Setup the halo box structure on this process.
     * (Private method since this is called as a helper method by the constructor.)
     *
     * Sets up the containers mHaloBoxes, mHalosRight, mHalosLeft
     */
    void SetupHaloBoxes();

    /** Needed for serialization **/
    friend class boost::serialization::access;

    /**
     * Serialize the box collection. It is possible to save and load on a different # of processes.
     *
     * @param archive the archive.
     * @param version the version number.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        //All methods over-ridden below.
    }

public:

    /**
     * @param boxWidth the width of each box (cut-off length in NodeBasedCellPopulation simulations)
     * @param domainSize the size of the domain, in the form (xmin, xmax, ymin, ymax) (etc)
     * @param isPeriodicInX whether the domain is periodic in the x direction
     * @param localRows the number of local rows in a parallel DistributedBoxCollection.
     *
     * Note that the domain size may be increased because each process should have at least one slice of boxes
     * in the largest dimension.  For example, if the box size is 1 and the domain is [(0,0,0), (3,3,3)] then,
     * if there are more than 3 processes the domain will be swollen to [(0,0,0), (3,3,num_procs)].  The
     * user is warned when this happens.
     */
    DistributedBoxCollection(double boxWidth, c_vector<double, 2*DIM> domainSize, bool isPeriodicInX = false, int localRows = PETSC_DECIDE);


    /**
     * Destructor - frees memory allocated to distributed vector.
     */
    ~DistributedBoxCollection();

    /**
     * Remove the list of nodes stored in each box.
     */
    void EmptyBoxes();

    /**
     * Update the halo boxes on this process, by transferring
     * the nodes to be sent into the lists mHaloNodesRight / Left.
     */
    void UpdateHaloBoxes();

    /**
     * @return the number of local rows / faces (2d, 3d) of boxes owned.
     */
    unsigned GetNumLocalRows() const;

    /**
     * @param globalIndex the global index of the box.
     * @return whether the box with global index globalIndex is owned by this process.
     */
    bool IsBoxOwned(unsigned globalIndex);

    /**
     * Get whether the box with global index globalIndex is interior on this process.
     * A box is interior if it doesn't share any boundary (even of zero length) with a halo box.
     *
     * @param globalIndex the global index of the box to check.
     * @return whether the box is interior or not.
     */
    bool IsInteriorBox(unsigned globalIndex);

    /**
     * @param globalIndex the global index of the box.
     * @return whether the box with global index globalIndex is a halo to this process.
     */
    bool IsHaloBox(unsigned globalIndex);

    /**
     * Given the (i,j,k) grid indices of a box in the collection, calculate its global index.
     *
     * @param gridIndices the (i,j,k) grid indices of the box
     * @return the global index of the box in the collection.
     */
    unsigned CalculateGlobalIndex(c_vector<unsigned, DIM> gridIndices);

    /**
     * @param pNode address of the node
     * @return the global index of the box that contains pNode.
     */
    unsigned CalculateContainingBox(Node<DIM>* pNode);

    /**
     * @param rLocation The point
     * @return the global index of the box that contains the rLocation.
     */
    unsigned CalculateContainingBox(c_vector<double, DIM>& rLocation);

    /**
     * Calculate x,y,z indices of box given its 'global' index.
     *
     * @param globalIndex the global index of the box
     * @return the grid indices (boxes across, boxes up, boxes deep) of a box
     */
    c_vector<unsigned, DIM> CalculateGridIndices(unsigned globalIndex);

    /**
     * If this box is out-of-bounds to the local process then it will attempt to return
     * a halo box (and trip an assertion is the box is completely out of scope.
     * @param boxIndex the index of the box to return
     * @return a reference to the box with global index boxIndex.
     */
    Box<DIM>& rGetBox(unsigned boxIndex);

    /**
     * Get a halo box.
     * @param boxIndex the index of the box to return
     * @return a reference to the halo box with global index boxIndex.
     */
    Box<DIM>& rGetHaloBox(unsigned boxIndex);

    /**
     * @return the total (global) number of boxes.
     */
    unsigned GetNumBoxes();

    /**
     * @return the number of locally owned boxes. Not including halos.
     */
    unsigned GetNumLocalBoxes();

    /**
     * @return #mDomainSize
     */
    c_vector<double, 2*DIM> rGetDomainSize() const;

    /**
     * @return Whether or not the local boxes have been set up.
     */
    bool GetAreLocalBoxesSet() const;

    /**
     * @return #mBoxWidth
     */
    double GetBoxWidth() const;

    /**
     * @return Whether the domain is periodic in x
     */
    bool GetIsPeriodicInX() const;

    /**
     * @return the number of rows in the DIM-1th direction on this process.
     */
    unsigned GetNumRowsOfBoxes() const;

    /**
     * A helper function to work out the optimal number of rows to be owned by this process, to balance the
     * number of nodes.
     *
     * @param localDistribution a vector containing the number of nodes in each row/face of boxes in 2d/3d
     * @return the updated number of rows, which will differ from current number by at most 2.
     */
    int LoadBalance(std::vector<int> localDistribution);

    /**
     *  Set up the local boxes (ie itself and its nearest-neighbours) for each of the boxes.
     *  This method just sets up half of the local boxes (for example, in 1D, local boxes for box0 = {1}
     *  local boxes for box1 = {2} not {0,2}, and so on. Similar to 2d, 3d.
     */
    void SetupLocalBoxesHalfOnly();

    /**
     * Set up the local boxes (ie itself and its nearest-neighbours) for each of the boxes.
     **/
    void SetupAllLocalBoxes();

    /**
     * Get the set of all the local boxes, i.e. itself and its nearest-neighbours.
     *
     * @param boxIndex the index of the box
     * @return the set containing the indices of boxes local to box boxIndex.  i.e. the box boxIndex itself and its nearest-neighbours.
     */
    std::set<unsigned>& rGetLocalBoxes(unsigned boxIndex);

    /**
     * @param pNode the node to test.
     * @return whether the point at pNode->rGetLocation() is owned on this process.
     */
    bool IsOwned(Node<DIM>* pNode);

    /**
     * @param location the location to test.
     * @return whether the point at location is owned on this process.
     */
    bool IsOwned(c_vector<double, DIM>& location);

    /**
     * Get the process that should own this node.
     * Currently only returns +/-1 of this process so assumes nodes don't move too far. //\ todo this should be fixed.
     *
     * @param pNode the node to be tested
     * @return the ID of the process that should own the node.
     */
    unsigned GetProcessOwningNode(Node<DIM>* pNode);

    /**
     * @return #mHaloNodesRight the list of nodes that are close to the right boundary
     */
    std::vector<unsigned>& rGetHaloNodesRight();

    /**
     * @return #mHaloNodesRight the list of nodes that are close to the left boundary
     */
    std::vector<unsigned>& rGetHaloNodesLeft();

    /**
     * Set whether to record node neighbour in the map rNodeNeighbours during CalculateNodePairs. Set to false for efficiency if not needed.
     *
     * @param calculateNodeNeighbours whether to store the neighbours.
     */
    void SetCalculateNodeNeighbours(bool calculateNodeNeighbours);

    /**
     *  Compute all the pairs of (potentially) connected nodes for cell_based simulations, ie nodes which are in a
     *  local box to the box containing the first node. **Note: the user still has to check that the node
     *  pairs are less than the cut-off distance apart.** The pairs are checked so that index1 < index2,
     *  so each connected pair of nodes is only in the set once.
     *
     *  @param rNodes all the nodes to be consider
     *  @param rNodePairs the return value, a set of pairs of nodes
     */
    void CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);

    /**
     *  The same as CalculateNodePairs() only we only work on boxes that are interior on this process. I.e. none of their local boxes are halo boxes.
     *
     *  @param rNodes all the nodes to be consider
     *  @param rNodePairs the return value, a set of pairs of nodes
     */
    void CalculateInteriorNodePairs(std::vector<Node<DIM>*>& rNodes, std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);

    /**
     *  The same as CalculateNodePairs() only we only work on boxes that are ''not'' interior on this process. I.e. some of their local boxes are halo boxes.
     *
     *  @param rNodes all the nodes to be consider
     *  @param rNodePairs the return value, a set of pairs of nodes
     */
    void CalculateBoundaryNodePairs(std::vector<Node<DIM>*>& rNodes, std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);

    /**
     * A method pulled out of CalculateNodePairs methods that adds all pairs of nodes from neighbouring boxes of the box with index boxIndex.
     *
     * @param boxIndex the box to add neighbours to.
     * @param rNodePairs the return value, a set of pairs of nodes
     */
    void AddPairsFromBox(unsigned boxIndex, std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs);

    /**
     * Calculate how many cells lie in each strip / face of boxes, used in load balancing
     *
     * @return A vector containing the number of nodes in each of the strips of boxes.
     */
    std::vector<int> CalculateNumberOfNodesInEachStrip();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DistributedBoxCollection)


namespace boost
{
namespace serialization
{
/**
 * Save information needed to reconstruct a box collection on load.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const DistributedBoxCollection<DIM> * t, const unsigned int file_version)
{
    // Save the number of rows that each process owns, so that on loading we can resume with
    // good load balance
    int num_local_rows = (int)(t->GetNumRowsOfBoxes());
    std::vector<int> num_rows(PetscTools::GetNumProcs());
    MPI_Gather(&num_local_rows, 1, MPI_INT, &num_rows[0], 1, MPI_INT, 0, PETSC_COMM_WORLD);

    if (PetscTools::AmMaster())
    {
        bool are_boxes_set = t->GetAreLocalBoxesSet();
        ar << are_boxes_set;

        c_vector<double, 2*DIM> domain_size = t->rGetDomainSize();
        for (unsigned i=0; i<2*DIM; i++)
        {
            ar << domain_size[i];
        }

        double box_width = t->GetBoxWidth();
        ar << box_width;

        unsigned num_procs = PetscTools::GetNumProcs();
        ar << num_procs;

        std::vector<int> const const_num_rows = num_rows;
        ar << const_num_rows;
    }
}
/**
 * De-serialize constructor parameters and initialize a DistributedBoxCollection.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, DistributedBoxCollection<DIM> * t, const unsigned int file_version)
{
    bool are_boxes_set;
    ar >> are_boxes_set;

    // Retrieve data from archive required to construct new instance of Node
    c_vector<double,2*DIM> domain_size;
    for (unsigned i=0; i<2*DIM; i++)
    {
        double coordinate;
        ar & coordinate; //resume coordinates one by one
        domain_size[i] = coordinate;
    }

    // Shrink the domain size slightly so we get the correct size on construction
    for (unsigned i=0; i<DIM; i++)
    {
        domain_size[2*i+1] -= 1e-14;
    }
    double cut_off;
    ar >> cut_off;

    unsigned num_original_procs;
    ar >> num_original_procs;

    int num_rows = PETSC_DECIDE;
    std::vector<int> original_rows;
    ar >> original_rows;
    if (num_original_procs == PetscTools::GetNumProcs())
    {
        num_rows = original_rows[PetscTools::GetMyRank()];
    }

    // Invoke inplace constructor to initialise instance. Assume non-periodic
    ::new(t)DistributedBoxCollection<DIM>(cut_off, domain_size, false, num_rows);

    if (are_boxes_set)
    {
        t->SetupLocalBoxesHalfOnly();
    }
}
}
} // namespace ...


#endif /*DISTRIBUTEDBOXCOLLECTION_HPP_*/
