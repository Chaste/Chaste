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

#ifndef MESHBASEDCELLPOPULATIONWITHGHOSTNODES_HPP_
#define MESHBASEDCELLPOPULATIONWITHGHOSTNODES_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "MeshBasedCellPopulation.hpp"

/**
 * A facade class encapsulating a mesh-based cell population with ghost nodes.
 *
 * If simulating a crypt with a mesh-based cell population, the mesh should be surrounded by at
 * least one layer of ghost nodes. These are nodes which do not correspond to a cell,
 * but are necessary for remeshing (because the remesher tries to create a convex hull
 * of the set of nodes) and visualization purposes. The MeshBasedCellPopulationWithGhostNodes
 * class deals with these ghost nodes, hiding the 'ghost nodes' concept from the
 * OffLatticeSimulation class, so the latter only ever deals with real cells.
 */
template<unsigned DIM>
class MeshBasedCellPopulationWithGhostNodes : public MeshBasedCellPopulation<DIM>
{
private:

    /** Just so that the test can test the private functions */
    friend class TestMeshBasedCellPopulationWithGhostNodes;

    /** Records whether a node is a ghost node or not */
    std::vector<bool> mIsGhostNode;

    /**
     * Spring stiffness for springs between ghost nodes.
     */
    double mGhostSpringStiffness;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This needs to be first so that MeshBasedCellPopulation::Validate() doesn't go mental.
        archive & mIsGhostNode;
        archive & mGhostSpringStiffness;
        archive & boost::serialization::base_object<MeshBasedCellPopulation<DIM, DIM> >(*this);
    }

    /**
     * Set the ghost nodes by taking in a set of which nodes indices are ghost nodes.
     *
     * @param rGhostNodeIndices set of node indices corresponding to ghost nodes
     */
    void SetGhostNodes(const std::set<unsigned>& rGhostNodeIndices);

    /**
     * This is called after a cell population has been constructed to check the
     * user gave consistent instructions. Check consistency of our
     * internal data structures:
     * Each node must have a cell associated with it OR must be a ghost node.
     *
     * It is called after cells are added or removed from MeshBasedCellPopulation
     * as it is an overridden virtual method.
     */
    void Validate();

protected:
    /**
     * Overridden method
     *
     * Calls #AcceptCellWriter across the whole population,
     * iterating in an appropriate way to skip ghost nodes.
     */
    virtual void AcceptCellWritersAcrossPopulation();

public:

    /**
     * Create a new cell population from a mesh and collection of cells.
     *
     * @param rMesh a mutable tetrahedral mesh
     * @param rCells cells corresponding to the nodes of the mesh
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     * @param ghostSpringStiffness spring stiffness used to move the ghost nodes defaults to 15.0.
     */
    MeshBasedCellPopulationWithGhostNodes(MutableMesh<DIM, DIM>& rMesh,
                                          std::vector<CellPtr>& rCells,
                                          const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                                          bool deleteMesh=false,
                                          double ghostSpringStiffness=15.0);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable tetrahedral mesh.
     * @param ghostSpringStiffness spring stiffness used to move the ghost nodes defaults to 15.0.
     */
    MeshBasedCellPopulationWithGhostNodes(MutableMesh<DIM, DIM>& rMesh,
                                          double ghostSpringStiffness=15.0);

    /**
     * Empty destructor so archiving works with static libraries.
     */
    virtual ~MeshBasedCellPopulationWithGhostNodes();

    /**
     * Overridden GetTetrahedralMeshForPdeModifier() method.
     *
     * @return a shared pointer to mrMesh.
     *
     * This method is called by AbstractGrowingDomainPdeModifier.
     */
    virtual TetrahedralMesh<DIM, DIM>* GetTetrahedralMeshForPdeModifier();

    /**
     * Overridden GetNeighbouringLocationIndices() method.
     *
     * Given a cell, returns the set of location indices corresponding to neighbouring cells.
     *
     * @param pCell a cell
     * @return the set of neighbouring location indices.
     */
    std::set<unsigned> GetNeighbouringLocationIndices(CellPtr pCell);

    /**
     * Applies the appropriate force to each ghost node in the population.
     * Called by AbstractNumericalMethod.
     */
    void ApplyGhostForces();

    /**
     * @return mIsGhostNode.
     */
    std::vector<bool>& rGetGhostNodes();

    /**
     * Overridden IsGhostNode() method.
     *
     * Find if a given node is a ghost node. The abstract method always returns false
     * but is overridden in subclasses.
     *
     * @param index the global index of a specified node
     *
     * @return whether the node is a ghost node
     */
    bool IsGhostNode(unsigned index);

    /**
     * @return the indices of those nodes that are ghost nodes.
     */
    std::set<unsigned> GetGhostNodeIndices();


    /**
     * Update mIsGhostNode if required by a remesh.
     *
     * @param rMap A map between node indices before and after remesh
     */
    void UpdateGhostNodesAfterReMesh(NodeMap& rMap);

    /**
     * This method is used to calculate the force between GHOST nodes.
     *
     * @param rNodeAGlobalIndex
     * @param rNodeBGlobalIndex
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, DIM> CalculateForceBetweenGhostNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex);

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the cell population and update mIsGhostNode.
     *
     * @param pNewCell  the cell to add
     * @param pParentCell pointer to a parent cell  - this is required for
     *  mesh-based cell populations
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    CellPtr AddCell(CellPtr pNewCell, CellPtr pParentCell);

    /**
     * Overridden OpenWritersFiles() method.
     *
     * Open all files in mCellPopulationWriters and mCellWriters for writing (not appending).
     *
     * @param rOutputFileHandler handler for the directory in which to open this file.
     */
    virtual void OpenWritersFiles(OutputFileHandler& rOutputFileHandler);

    /**
     * Overridden WriteVtkResultsToFile() method.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     */
    virtual void WriteVtkResultsToFile(const std::string& rDirectory);

    /**
     * Outputs CellPopulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile);
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedCellPopulationWithGhostNodes)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MeshBasedCellPopulationWithGhostNodes.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MeshBasedCellPopulationWithGhostNodes<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const MutableMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a MeshBasedCellPopulationWithGhostNodes.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MeshBasedCellPopulationWithGhostNodes<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    MutableMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)MeshBasedCellPopulationWithGhostNodes<DIM>(*p_mesh);
}
}
} // namespace

#endif /*MESHBASEDCELLPOPULATIONWITHGHOSTNODES_HPP_*/
