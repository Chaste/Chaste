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

#ifndef NODEBASEDCELLPOPULATION_HPP_
#define NODEBASEDCELLPOPULATION_HPP_

#include "AbstractCentreBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "BoxCollection.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A NodeBasedCellPopulation is a CellPopulation consisting of only nodes in space with associated cells.
 * There are no elements and no mesh.
 */
template<unsigned DIM>
class NodeBasedCellPopulation : public AbstractCentreBasedCellPopulation<DIM>
{
    friend class TestNodeBasedCellPopulation;
    friend class TestBoxCollection;

protected:

    /** Reference to the mesh. */
    NodesOnlyMesh<DIM>& mrMesh;

private:

    /** Pointer to a Node box collection. */
    BoxCollection<DIM>* mpBoxCollection;

    /** Vector of minimal spatial positions in each dimension. */
    c_vector<double, DIM> mMinSpatialPositions;

    /** Vector of maximal spatial positions in each dimension. */
    c_vector<double, DIM> mMaxSpatialPositions;

    /** Node pairs for force calculations. */
    std::set< std::pair<Node<DIM>*, Node<DIM>* > > mNodePairs;

    /** Whether to delete the nodes-only mesh (taken in one of the constructors, defaults to false). */
    bool mDeleteMesh;

    /**
     * Mechanics cut off length.
     * Used in order to calculate the BoxCollection.
     */
    double mMechanicsCutOffLength;

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
        archive & mMechanicsCutOffLength;

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
    /**
     * Move the node with a given index to a new point in space.
     *
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation);

private:
    /**
     * Check consistency of our internal data structures.
     */
    void Validate();

    /**
     * Method for Initially Splitting up cell population into neighbouring boxes, to decrease runtime.
     *
     * @param cutOffLength length of spring cut off between nodes
     * @param domainSize c_vector of size 2*dimension reads minX, maxX, minY, maxY, etc
     */
    void SplitUpIntoBoxes(double cutOffLength, c_vector<double, 2*DIM> domainSize);

    /**
     * Loops over nodes and sets mMinSpatialPositions and mMaxSpatialPositions
     */
    void FindMaxAndMin();

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
     * Default constructor.
     *
     * Note that the cell population will take responsibility for freeing the memory used by the nodes.
     *
     * @param rMesh a mutable nodes-only mesh
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh whether to delete nodes-only mesh in destructor
     */
    NodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh,
                            std::vector<CellPtr>& rCells,
                            const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                            bool deleteMesh=false);

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
    ~NodeBasedCellPopulation();

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
     * Reset the member variables #mNodePairs and #mpBoxCollection.
     */
    void Clear();

    /**
     * Remove nodes that have been marked as deleted and update the node cell map.
     *
     * @param hasHadBirthsOrDeaths whether cell population has had Births Or Deaths
     */
    void Update(bool hasHadBirthsOrDeaths=true);

    /**
     * @return pointer to a node box collection.
     */
    BoxCollection<DIM>* GetBoxCollection();

    /**
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
     * @return mMechanicsCutOffLength
     */
    double GetMechanicsCutOffLength();

    /**
     * Set mMechanicsCutOffLength.
     *
     * @param mechanicsCutOffLength  the new value of mMechanicsCutOffLength
     */
    void SetMechanicsCutOffLength(double mechanicsCutOffLength);

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
     * Overridden GetNeighbouringNodeIndices() method.
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
     *  mesh-based cell populations
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    virtual CellPtr AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell=CellPtr());

    /**
     * Overridden GetVolumeOfCell() method.
     * 
     * @param pCell boost shared pointer to a cell
     */
    double GetVolumeOfCell(CellPtr pCell);
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
