/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef VERTEXBASEDTISSUE_HPP_
#define VERTEXBASEDTISSUE_HPP_

#include "AbstractTissue.hpp"
#include "VertexMesh.hpp"
#include "ArchiveLocationInfo.hpp"

#include <climits> // work around boost bug

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>


/**
 * A facade class encapsulating a vertex-based 'tissue'.
 *
 * Contains a group of cells and maintains the associations
 * between TissueCells and elements in the VertexMesh.
 *
 */
template<unsigned DIM>
class VertexBasedTissue : public AbstractTissue<DIM>
{
private:

    /** Vertex-based mesh associated with the tissue. */
    VertexMesh<DIM, DIM>& mrMesh;

    /** A cache of where the results are going (used for VTK writer). */
    std::string mDirPath;

    /** Meta results file for VTK. */
    out_stream mpVtkMetaFile;

    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this tissue has been de-serialized.
     */
    bool mDeleteMesh;

    /** Results file for elements. */
    out_stream mpElementFile;

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
        archive & boost::serialization::base_object<AbstractTissue<DIM> >(*this);
    }

    /**
     * Check the consistency of internal data structures.
     * Each VertexElement must have a TissueCell associated with it.
     */
    void Validate();

public:

    /**
     * Create a new tissue facade from a mesh and collection of cells.
     *
     * There must be precisely one TissueCell for each VertexElement in
     * the mesh.
     *
     * @param rMesh reference to a VertexMesh
     * @param rCells reference to a vector of TissueCells
     * @param deleteMesh set to true if you want the tissue to free the mesh memory on destruction
     * @param validate whether to validate the tissue when it is created (defaults to true)
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    VertexBasedTissue(VertexMesh<DIM, DIM>& rMesh,
                      const std::vector<TissueCell>& rCells,
                      bool deleteMesh=false,
                      bool validate=true,
                      const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a vertex mesh.
     */
    VertexBasedTissue(VertexMesh<DIM, DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~VertexBasedTissue();

    /**
     * Overridden GetDampingConstant() method.
     *
     * @param nodeIndex the global index of this node
     * @return the average damping constant of the cells surrounding the node.
     */
    double GetDampingConstant(unsigned nodeIndex);

    /**
     * Get the adhesion parameter for the edge between two given nodes.
     *
     * \todo This method should be changed/overridden if we require differential adhesion
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     *
     * @return the adhesion parameter for this edge.
     */
    double GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB);

    /**
     * @return reference to  mrMesh.
     */
    VertexMesh<DIM, DIM>& rGetMesh();

    /**
     * @return const reference to mrMesh (used in archiving).
     */
    const VertexMesh<DIM, DIM>& rGetMesh() const;

    /**
     * Get a particular VertexElement.
     *
     * @param elementIndex the global index of the VertexElement
     *
     * @return a pointer to the VertexElement.
     */
    VertexElement<DIM, DIM>* GetElement(unsigned elementIndex);

    /**
     * @return the number of VertexElements in the tissue.
     */
    unsigned GetNumElements();

    /**
     * Overridden GetNumNodes() method.
     *
     * @return the number of nodes in the tissue.
     */
    unsigned GetNumNodes();

    /**
     * Overridden GetLocationOfCellCentre() method.
     * Find where a given cell is in space.
     *
     * \todo If required, we could come up with a more clever definition of cell location
     *       for a VertexTissue (for example, there is no guarantee of convexity so the
     *       centre of mass may lie outside the element)
     *
     * @param rCell the cell
     *
     * @return the location of the centre of mass of the element corresponding to this cell.
     */
    c_vector<double, DIM> GetLocationOfCellCentre(TissueCell& rCell);

    /**
     * Overridden GetNode() method.
     *
     * @param index global index of the specified node
     *
     * @return a pointer to the node.
     */
    Node<DIM>* GetNode(unsigned index);

    /**
     * Overridden AddNode() method.
     *
     * Add a new node to the tissue.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in tissue
     */
    unsigned AddNode(Node<DIM>* pNewNode);

    /**
     * Overridden UpdateNodeLocations() method.
     *
     * @param rNodeForces a vector containing the force on each node in the tissue
     * @param dt the time step
     */
    void UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt);

    /**
     * Overridden SetNode() method.
     *
     * Move the node with a given index to a new point in space.
     *
     * @param index the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned index, ChastePoint<DIM>& rNewLocation);

    /**
     * Get a pointer to the element corresponding to a given TissueCell.
     *
     * @param rCell the cell
     *
     * @return pointer to the element.
     */
    VertexElement<DIM, DIM>* GetElementCorrespondingToCell(TissueCell& rCell);

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the tissue.
     *
     * @param rNewCell  the cell to add
     * @param rCellDivisionVector  if this vector has any non-zero component, then it is used as the axis
     *     along which the parent cell divides
     * @param pParentCell pointer to a parent cell (if required)
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell* AddCell(TissueCell& rNewCell, const c_vector<double,DIM>& rCellDivisionVector, TissueCell* pParentCell=NULL);

    /**
     * Remove all cells labelled as dead.
     *
     * Note that after calling this method the tissue will be in an inconsistent state until
     * the equivalent of a 'remesh' is performed! So don't try iterating over cells or anything
     * like that.
     *
     * @return number of cells removed
     */
    unsigned RemoveDeadCells();

    /**
     * Overridden IsCellAssociatedWithADeletedLocation() method.
     *
     * @param rCell the cell
     * @return whether a given cell is associated with a deleted element.
     */
    bool IsCellAssociatedWithADeletedLocation(TissueCell& rCell);

    /**
     * Remove the VertexElements which have been marked as deleted, perform
     * any cell rearrangements if required, and update the correspondence
     * with TissueCells.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether tissue has had Births Or Deaths
     * not needed in this tissue class
     */
    void Update(bool hasHadBirthsOrDeaths=true);

    /**
     * Get the target area of a given cell. This grows linearly from
     * 0.5*A to A during the G1 phase of the cell cycle, then remains
     * at A for the rest of the cell cycle, where A denotes the TissueConfig
     * member variable mMatureCellTargetArea.
     *
     * @param rCell the cell
     * @return the cell's target area
     */
    double GetTargetAreaOfCell(const TissueCell& rCell);

    /**
     * Overridden CreateOutputFiles() method.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     * @param cleanOutputDirectory  whether to delete the contents of the output directory prior to output file creation
     */
    void CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory);
    /**
     * Overridden CloseOutputFiles() method.
     */
    void CloseOutputFiles();

    /**
     * Overridden WriteResultsToFiles() method.
     */
    void WriteResultsToFiles();

    /**
     * Overridden GenerateCellResultsAndWriteToFiles() method.
     */
    virtual void GenerateCellResultsAndWriteToFiles();

};

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedTissue)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VertexBasedTissue.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VertexBasedTissue<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const VertexMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a VertexBasedTissue.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VertexBasedTissue<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    VertexMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)VertexBasedTissue<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*VERTEXBASEDTISSUE_HPP_*/

