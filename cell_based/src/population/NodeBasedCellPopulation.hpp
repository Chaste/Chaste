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

#ifndef NODEBASEDCELLPOPULATION_HPP_
#define NODEBASEDCELLPOPULATION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCentreBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "BoxCollection.hpp"

/**
 * A NodeBasedCellPopulation is a CellPopulation consisting of only nodes in space with associated cells.
 * There are no elements and no mesh.
 */
template<unsigned DIM>
class NodeBasedCellPopulation : public AbstractCentreBasedCellPopulation<DIM>
{
    friend class TestNodeBasedCellPopulation;
    friend class TestNodeBasedCellPopulationParallelMethods;
    friend class TestBoxCollection;

private:

    /** Static cast of the mesh from AbstractCellPopulation. */
    NodesOnlyMesh<DIM>* mpNodesOnlyMesh;

    /** Vector of minimal spatial positions in each dimension. */
    c_vector<double, DIM> mMinSpatialPositions;

    /** Vector of maximal spatial positions in each dimension. */
    c_vector<double, DIM> mMaxSpatialPositions;

    /** Node pairs for force calculations. */
    std::set< std::pair<Node<DIM>*, Node<DIM>* > > mNodePairs;

    /** Indices of neighbours of local nodes. */
    std::map<unsigned, std::set<unsigned> > mNodeNeighbours;

    /** Whether to delete the nodes-only mesh (taken in one of the constructors, defaults to false). */
    bool mDeleteMesh;

    /**
     * Whether or not to have cell radii updated from CellData defaults to false.
     */
    bool mUseVariableRadii;

    /** The cells to send to the right process */
    std::set<std::pair<CellPtr, Node<DIM>*> > mCellsToSendRight;

    /** The cells to send to the left process */
    std::set<std::pair<CellPtr, Node<DIM>*> > mCellsToSendLeft;

    /** A pointer to the cells received from the right process */
    std::set<std::pair<CellPtr, Node<DIM>*> >* mpCellsRecvRight;

    /** A pointer to the cells received from the left process */
    std::set<std::pair<CellPtr, Node<DIM>*> >* mpCellsRecvLeft;

    /** The tag used to send and recieve cell information */
    static const unsigned mCellCommunicationTag = 123;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the nodes is handled by load/save_construct_data,
     * so we don't actually have to do anything here except delegate to the base class.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCentreBasedCellPopulation<DIM> >(*this);
        archive & mUseVariableRadii;

        this->Validate();
    }

    /**
     * Overridden AddNode() method.
     *
     * Add a new node to the cell population.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in cell population
     */
    unsigned AddNode(Node<DIM>* pNewNode);


protected:

#undef COVERAGE_IGNORE // Avoid prototypes being treated as code by gcov
    /**
     * Update mIsParticle if required by a remesh.
     *
     * @param rMap A map between node indices before and after remesh
     */
    virtual void UpdateParticlesAfterReMesh(NodeMap& rMap);

    /**
     * Check consistency of our internal data structures.
     */
    virtual void Validate();

private:

    /**
     * Overridden WriteCellVolumeResultsToFile() method.
     */
    void WriteCellVolumeResultsToFile();

    /**
     * Overridden WriteVtkResultsToFile() method.
     */
    void WriteVtkResultsToFile();

public:

    /**
     * Move the node with a given index to a new point in space.
     *
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation);

    /**
     * Default constructor.
     *
     * Note that the cell population will take responsibility for freeing the memory used by the nodes.
     *
     * @param rMesh a mutable nodes-only mesh
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh whether to delete nodes-only mesh in destructor
     * @param validate whether to call Validate() in the constructor or not
     */
    NodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh,
                            std::vector<CellPtr>& rCells,
                            const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                            bool deleteMesh=false,
                            bool validate=true);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable nodes-only mesh
     */
    NodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh);

    /**
     * Destructor.
     *
     * Frees all our node memory.
     */
    virtual ~NodeBasedCellPopulation();

    /**
     * @return reference to  mrMesh.
     */
    NodesOnlyMesh<DIM>& rGetMesh();

    /**
     * @return const reference to mrMesh (used in archiving).
     */
    const NodesOnlyMesh<DIM>& rGetMesh() const;

    /**
     * @return the number of nodes in the cell population.
     */
    unsigned GetNumNodes();

    /**
     * Overridden GetNode() method.
     *
     * @param index global index of the specified node
     *
     * @return a pointer to the node with a given index.
     */
    Node<DIM>* GetNode(unsigned index);

    /**
     * Remove all cells labelled as dead.
     *
     * Note that after calling this method the cell population will be in an inconsistent state until
     * the equivalent of a 'remesh' is performed! So don't try iterating over cells or anything
     * like that.
     *
     * @return number of cells removed
     */
    unsigned RemoveDeadCells();

    /**
     * Reset the member variables #mNodePairs and mpBoxCollection in the underlying mesh.
     */
    void Clear();

    /**
     * Remove nodes that have been marked as deleted and update the node cell map.
     *
     * @param hasHadBirthsOrDeaths whether cell population has had Births Or Deaths
     */
    void Update(bool hasHadBirthsOrDeaths=true);

    /**
     * Overridden rGetNodePairs method
     *
     * @return Node pairs for force calculation.
     */
    std::set< std::pair<Node<DIM>*, Node<DIM>* > >& rGetNodePairs();

    /**
     * Outputs CellPopulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile);

    /**
     * @return the maximum interaction distance between cells, defined in NodesOnlyMesh.
     */
    double GetMechanicsCutOffLength();

    /**
     * @return mUseVariableRadii
     */
    bool GetUseVariableRadii();

    /**
     * Set mUseVariableRadii.
     *
     * @param useVariableRadii the new value of mUseVariableRadii
     */
    void SetUseVariableRadii(bool useVariableRadii=true);

    /**
     * Overridden GetWidth() method.
     *
     * Calculate the 'width' of any dimension of the cell population by computing
     * the maximum distance between any nodes in this dimension.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    double GetWidth(const unsigned& rDimension);

    /**
     * Overridden GetSizeOfCellPopulation to work in parallel.
     *
     * @return the size of the cell population.
     */
    c_vector<double, DIM> GetSizeOfCellPopulation();

    /**
     * Overridden GetNeighbouringNodeIndices() method.
     *
     * Only returns nodes that are strictly touching each other.
     *
     * @param index the node index
     * @return the set of neighbouring node indices.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index);

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the cell population and update the vector of cell radii in the NodesOnlyMesh.
     *
     * @param pNewCell  the cell to add
     * @param rCellDivisionVector  the position in space at which to put it
     * @param pParentCell pointer to a parent cell - this is required for
     *  node-based cell populations
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    virtual CellPtr AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell);

    /**
     * Overridden GetVolumeOfCell() method.
     *
     * @param pCell boost shared pointer to a cell
     * @return volume via associated mesh node
     */
    double GetVolumeOfCell(CellPtr pCell);

    /////////////////////////////////////////////////////
    // Parallel methods
    /////////////////////////////////////////////////////

    /**
     * Send the contents of #mCellsToSendRight/Left to
     * neighbouring processes and receive from them into
     * #mpCellsRecvRight/Left.
     */
    void SendCellsToNeighbourProcesses();

    /**
     * Helper method to find and pack up nodes and cells together
     *
     * @param nodeIndex the global index of the node.
     * @return the pair.
     */
    std::pair<CellPtr, Node<DIM>*> GetCellNodePair(unsigned nodeIndex);

    /**
     * Add the node and cell with index nodeIndex to the list of cells to send
     * to the process right.
     *
     * @param nodeIndex the index of the node and cell to send.
     */
    void AddNodeAndCellToSendRight(unsigned nodeIndex);

    /**
     * Add the node and cell with index nodeIndex to the list of cells to send
     * to the process left.
     *
     * @param nodeIndex the index of the node and cell to send.
     */
    void AddNodeAndCellToSendLeft(unsigned nodeIndex);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a NodeBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const NodeBasedCellPopulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const NodesOnlyMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a NodeBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, NodeBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    NodesOnlyMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)NodeBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*NODEBASEDCELLPOPULATION_HPP_*/
