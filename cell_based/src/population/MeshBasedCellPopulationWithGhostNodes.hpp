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

#ifndef MESHBASEDCELLPOPULATIONWITHGHOSTNODES_HPP_
#define MESHBASEDCELLPOPULATIONWITHGHOSTNODES_HPP_

#include "MeshBasedCellPopulation.hpp"
#include "TrianglesMeshReader.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

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
        archive & boost::serialization::base_object<MeshBasedCellPopulation<DIM> >(*this);
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
     * Overridden UpdateNodeLocation() method.
     *
     * Update the location of each node in the cell population given
     * a vector of forces on nodes and a time step over which
     * to integrate the equations of motion.
     *
     * @param rNodeForces  forces on nodes
     * @param dt  time step
     */
    void UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt);

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
     * Update the GhostNode positions using the spring force model with rest length=1.
     * Forces are applied to ghost nodes from connected ghost and normal nodes.
     *
     * @param dt
     */
    void UpdateGhostPositions(double dt);

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
     * @param rCellDivisionVector  the position in space at which to put it
     * @param pParentCell pointer to a parent cell  - this is required for
     *  mesh-based cell populations
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    CellPtr AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell);

    /**
     * Overridden WriteVtkResultsToFile method.
     */
    void WriteVtkResultsToFile();

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
    Archive & ar, const MeshBasedCellPopulationWithGhostNodes<DIM> * t, const BOOST_PFTO unsigned int file_version)
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
