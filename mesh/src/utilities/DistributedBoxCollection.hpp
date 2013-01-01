/*

Copyright (c) 2005-2013, University of Oxford.
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

    /** A set of nodes that are on neighbouring processes */
    std::set<Node<DIM>* > mHaloNodes;

    /** Set of Nodes that are halos of adjacent right process, but lie locally */
    std::set<unsigned> mHaloNodesRight;

    /** Set of Nodes that are halos of adjacent left process, but lie locally */
    std::set<unsigned> mHaloNodesLeft;

    /** Map of global to local indices of boxes. **/
    std::map<unsigned, unsigned> mBoxesMapping;

    /** Map of global to local indices of halo boxes in mHaloBoxes. **/
    std::map<unsigned, unsigned> mHaloBoxesMapping;

    /** The domain being partitioned. */
    c_vector<double, 2*DIM> mDomainSize;

    /** The local portion of the domain. */
    c_vector<double, 2*DIM> mLocalDomainSize;

    /** The width of each box. */
    double mBoxWidth;

    /** The total number of boxes across all processes */
    unsigned mNumBoxes;

    /** Number of boxes in each direction. */
    c_vector<unsigned, DIM> mNumBoxesEachDirection;

    /** The boxes local (itself and nearest neighbour) to a given box. */
    std::vector< std::set<unsigned> > mLocalBoxes;

    /** The smallest index of the boxes owned by this process. */
    unsigned mMinBoxIndex;

    /** The largest index of the boxes owned by this process. */
    unsigned mMaxBoxIndex;

    /** A fudge (box swelling) factor to deal with 32-bit floating point issues. */
    const static double mFudge = 5e-14;

    /** Whether the domain is periodic in the X dimension Note this currently only works for DIM=2.*/
    bool mIsPeriodicInX;

    /** Whether the local boxes have been setup or not. */
    bool mAreLocalBoxesSet;

    /** The rank of the process to the right */
    unsigned mProcRight;

    /** The rank of the process to the right */
    unsigned mProcLeft;

    /** A distributed vector that govens ownership of rows of boxes */
    DistributedVector* mpDistributedBoxStacks;

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
     * Constructor
     *
     * @param boxWidth the width of each box (cut-off length in NodeBasedCellPopulation simulations)
     * @param domainSize the size of the domain, in the form (xmin, xmax, ymin, ymax) (etc)
     * @param isPeriodicInX whether the domain is peiodic in the x direction
     * @param localRows the number of local rows in a parallel DistributedBoxCollection.
     */
    DistributedBoxCollection(double boxWidth, c_vector<double, 2*DIM> domainSize, bool isPeriodicInX = false, int localRows = PETSC_DECIDE);


    /**
     * Destructor
     */
    ~DistributedBoxCollection();

    /**
     * Setup the halo box structure on this process.
     *
     * Sets up the containers mHaloBoxes, mHalosRight, mHalosLeft
     */
    void SetupHaloBoxes();

    /**
     * Update the halo boxes on this processes, by transferring
     * the nodes to be sent into the lists mHaloNodesRight / Left
     */
    void UpdateHaloBoxes();

    /**
     * Get whether the box with global index globalIndex is owned by this process.
     *
     * @param globalIndex the global index of the box.
     * @return whether the box is owned.
     */
    bool GetBoxOwnership(unsigned globalIndex);

    /**
     * Get whether the box with global index globalIndex is a halo to this process.
     *
     * @param globalIndex the global index of the box.
     * @return whether the box is halo-owned.
     */
    bool GetHaloBoxOwnership(unsigned globalIndex);

    /**
     * Given the (i,j,k) co-ordniates of a box in the collection, calculate its global index.
     *
     * @param coordinateIndices the co-ordinate indices of the box (across, up, deep)
     * @return the global index of the box in the collection.
     */
    unsigned CalculateGlobalIndex(c_vector<unsigned, DIM> coordinateIndices);

    /**
     * Calculate which box this node is contained in.
     *
     * @param pNode address of the node
     * @return the global index of the box that contains pNode.
     */
    unsigned CalculateContainingBox(Node<DIM>* pNode);

    /**
     * Calculate which box a point is contained in
     *
     * @param rLocation The point
     * @return the global index of the box that contains the location rLocation.
     */
    unsigned CalculateContainingBox(c_vector<double, DIM>& rLocation);

    /**
     * Calculate x,y,z indices of box given its 'global' index.
     *
     * @param globalIndex the global index of the box
     * @return the co-ordinate indicies (boxes across, boxes up, boxes deep) of a box
     */
    c_vector<unsigned, DIM> CalculateCoordinateIndices(unsigned globalIndex);

    /**
     * Get a box.
     * @param boxIndex the index of the box to return
     * @return a reference to the box.
     */
    Box<DIM>& rGetBox(unsigned boxIndex);

    /**
     * Get the number of boxes.
     *
     * @return mNumBoxes
     */
    unsigned GetNumBoxes();

    /**
     * Get the number of locally owned boxes.
     *
     * @return the number of locally owned boxes. Not including HaloBoxes.
     */
    unsigned GetNumLocalBoxes();

    /**
     * Get the global domain size. Used in serialisation.
     *
     * @return mDomainSize
     */
    c_vector<double, 2*DIM> GetDomainSize() const;

    /**
     * Return whether the local boxes have been set up or not. Used in serialisation.
     *
     * @return mAreLocalBoxesSet
     */
    bool GetAreLocalBoxesSet() const;

    /**
     * Return the width of each box. Used in serialisation.
     *
     * @return mBoxWdith
     */
    double GetBoxWidth() const;

    /**
     * Return the number of rows in the DIM-1th direction
     * Used in serialization
     *
     * @return the number of rows in the DIM-1th direction on this process.
     */
    unsigned GetNumRowsOfBoxes() const;

    /**
     *  Set up the local boxes (ie itself and its nearest-neighbours) for each of the boxes.
     *  Just set up half of the local boxes (for example, in 1D, local boxes for box0 = {1}
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
     * @return the set containing the indices of boxes local to box boxIndex.
     */
    std::set<unsigned> GetLocalBoxes(unsigned boxIndex);

    /**
     * Return the local portion of space owned by this process
     *
     * @return mLocalDomainSize
     */
    c_vector<double, 2*DIM> GetLocalDomainSize();

    /**
     *  Compute all the pairs of (potentially) connected nodes for cell_based simulations, ie nodes which are in a
     *  local box to the box containing the first node. **Note: the user still has to check that the node
     *  pairs are less than the cut-off distance apart.** The pairs are checked so that index1 < index2,
     *  so each connected pair of nodes is only in the set once.
     *
     *  @param rNodes all the nodes to be consider
     *  @param rNodePairs the return value, a set of pairs of nodes
     *  @param rNodeNeighbours the other return value, the neighbours of each node.
     */
    void CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::set<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs, std::map<unsigned, std::set<unsigned> >& rNodeNeighbours);
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
    std::vector<int> num_rows;
    if (PetscTools::AmMaster())
    {
        num_rows.resize(PetscTools::GetNumProcs());
    }
    MPI_Gather(&num_local_rows, 1, MPI_INT, &num_rows[0], 1, MPI_INT, 0, PETSC_COMM_WORLD);
    
    if (PetscTools::AmMaster())
    {
        bool are_boxes_set = t->GetAreLocalBoxesSet();
        ar << are_boxes_set;

        c_vector<double, 2*DIM> domain_size = t->GetDomainSize();

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
        t->SetupHaloBoxes();
        t->SetupLocalBoxesHalfOnly();
    }
}
}
} // namespace ...


#endif /*DISTRIBUTEDBOXCOLLECTION_HPP_*/
