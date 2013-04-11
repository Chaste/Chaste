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

#ifndef VERTEXBASEDCELLPOPULATION_HPP_
#define VERTEXBASEDCELLPOPULATION_HPP_

#include "AbstractOffLatticeCellPopulation.hpp"
#include "MutableVertexMesh.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"

/**
 * A facade class encapsulating a vertex-based cell population.
 *
 * Contains a group of cells and maintains the associations
 * between CellPtrs and elements in the MutableVertexMesh.
 *
 */
template<unsigned DIM>
class VertexBasedCellPopulation : public AbstractOffLatticeCellPopulation<DIM>
{
private:

    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this cell population has been de-serialized.
     */
    bool mDeleteMesh;

    /**
     * A static cast of the AbstractMesh from AbstractCellPopulation
     * for use in this class
     */
    MutableVertexMesh<DIM, DIM>* mpMutableVertexMesh;

    /** Results file for elements. */
    out_stream mpVizElementsFile;

    /** Results file for locations of T1Swaps. */
    out_stream mpT1SwapLocationsFile;

    /** Results file for locations of T3Swaps. */
    out_stream mpT3SwapLocationsFile;

    /** Whether to output the locations of T1Swaps and T3Swaps to files. Defaults to true. */
    bool mOutputCellRearrangementLocations;

    /**
     * Overridden WriteVtkResultsToFile() method.
     */
    void WriteVtkResultsToFile();

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
        archive & boost::serialization::base_object<AbstractOffLatticeCellPopulation<DIM> >(*this);
        archive & mOutputCellRearrangementLocations;
    }

    /**
     * Check the consistency of internal data structures.
     * Each VertexElement must have a CellPtr associated with it.
     */
    void Validate();

public:

    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely one CellPtr for each VertexElement in
     * the mesh.
     *
     * @param rMesh reference to a
     * @param rCells reference to a vector of CellPtrs
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     * @param validate whether to validate the cell population when it is created (defaults to true)
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    VertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh,
                              std::vector<CellPtr>& rCells,
                              bool deleteMesh=false,
                              bool validate=true,
                              const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a vertex mesh.
     */
    VertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~VertexBasedCellPopulation();

    /**
     * Overridden GetDampingConstant() method.
     *
     * @param nodeIndex the global index of this node
     * @return the average damping constant of the cells surrounding the node.
     */
    double GetDampingConstant(unsigned nodeIndex);

    /**
     * @return reference to  mrMesh.
     */
    MutableVertexMesh<DIM, DIM>& rGetMesh();

    /**
     * @return const reference to mrMesh (used in archiving).
     */
    const MutableVertexMesh<DIM, DIM>& rGetMesh() const;

    /**
     * Get a particular VertexElement.
     *
     * @param elementIndex the global index of the VertexElement
     *
     * @return a pointer to the VertexElement.
     */
    VertexElement<DIM, DIM>* GetElement(unsigned elementIndex);

    /**
     * @return the number of VertexElements in the cell population.
     */
    unsigned GetNumElements();

    /**
     * Overridden GetNumNodes() method.
     *
     * @return the number of nodes in the cell population.
     */
    unsigned GetNumNodes();

    /**
     * Overridden GetLocationOfCellCentre() method.
     *
     * Find the centre of mass of a given cell (assuming uniform density).
     * Note that, as there is no guarantee of convexity, this may lie
     * outside the VertexElement corresponding to the cell.
     *
     * @param pCell a cell in the population
     *
     * @return the location of the centre of mass of the element corresponding to this cell.
     */
    c_vector<double, DIM> GetLocationOfCellCentre(CellPtr pCell);

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
     * Add a new node to the cell population.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in cell population
     */
    unsigned AddNode(Node<DIM>* pNewNode);

    /**
     * Overridden UpdateNodeLocations() method.
     *
     * @param dt the time step
     */
    void UpdateNodeLocations(double dt);

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
     * Get a pointer to the element corresponding to a given CellPtr.
     *
     * @param pCell the cell
     *
     * @return pointer to the element.
     */
    VertexElement<DIM, DIM>* GetElementCorrespondingToCell(CellPtr pCell);

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the cell population.
     *
     * @param pNewCell  the cell to add
     * @param rCellDivisionVector  if this vector has any non-zero component, then it is used as the axis
     *     along which the parent cell divides
     * @param pParentCell pointer to a parent cell (if required)
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    CellPtr AddCell(CellPtr pNewCell,
                    const c_vector<double,DIM>& rCellDivisionVector,
                    CellPtr pParentCell=CellPtr());

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
     * Overridden IsCellAssociatedWithADeletedLocation() method.
     *
     * @param pCell the cell
     * @return whether a given cell is associated with a deleted element.
     */
    bool IsCellAssociatedWithADeletedLocation(CellPtr pCell);

    /**
     * Remove the VertexElements which have been marked as deleted, perform
     * any cell rearrangements if required, and update the correspondence
     * with CellPtrs.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether cell population has had Births Or Deaths
     * not needed in this cell population class
     */
    void Update(bool hasHadBirthsOrDeaths=true);

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
     * Write the current index and location (of the centre) of each element in mrMesh, as well as the ID and
     * the area (in 2D) or volume (in 3D) of its corresponding cell, to mpCellVolumesFile.
     */
    void WriteCellVolumeResultsToFile();

    /**
     * Overridden GetVolumeOfCell() method.
     *
     * @param pCell boost shared pointer to a cell
     * @return volume via associated mesh element
     */
    double GetVolumeOfCell(CellPtr pCell);

    /**
     * Overridden GenerateCellResultsAndWriteToFiles() method.
     */
    virtual void GenerateCellResultsAndWriteToFiles();

    /**
     * @return mOutputCellRearrangementLocations
     */
    bool GetOutputCellRearrangementLocations();

    /**
     * Set mOutputCellRearrangementLocations.
     *
     * @param outputCellRearrangementLocations the new value of mOutputCellRearrangementLocations
     */
    void SetOutputCellRearrangementLocations(bool outputCellRearrangementLocations);

    /**
     * Overridden OutputCellPopulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile);

    /**
     * Overridden GetWidth() method.
     *
     * Calculate the 'width' of any dimension of the cell population by calling
     * GetWidth() on the mesh.
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
     * @return a tetrahedral mesh using the nodes in the VertexMesh
     * as well as an additional node at the centre of each VertexElement.
     * This method only works in 2D.
     */
    TetrahedralMesh<DIM, DIM>* GetTetrahedralMeshUsingVertexMesh();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VertexBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VertexBasedCellPopulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const MutableVertexMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a VertexBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VertexBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    MutableVertexMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)VertexBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*VERTEXBASEDCELLPOPULATION_HPP_*/

