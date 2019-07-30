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

#ifndef POTTSBASEDCELLPOPULATION_HPP_
#define POTTSBASEDCELLPOPULATION_HPP_

#include "AbstractOnLatticeCellPopulation.hpp"
#include "PottsMesh.hpp"
#include "VertexMesh.hpp"
#include "AbstractUpdateRule.hpp"
#include "MutableMesh.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A facade class encapsulating a cell population under the Cellular
 * Potts Model framework.
 *
 * Contains a group of cells and maintains the associations
 * between CellPtrs and elements in a specialised PottsMesh class.
 *
 * The code currently requires the PottsMesh object to be fixed,
 * in the sense that no new nodes or elements can be added.
 */
template<unsigned DIM>
class PottsBasedCellPopulation : public AbstractOnLatticeCellPopulation<DIM>
{
    friend class TestPottsBasedCellPopulation;

private:

    /**
     * Pointer to a VertexMesh object that stores the Element tessellation that is used to
     * visualize mrMesh. The tessellation is created by calling CreateElementTessellation()
     * and can be accessed by calling GetElementTessellation().
     */
    VertexMesh<DIM,DIM>* mpElementTessellation;

    /**
     * A static cast of the Abstract mesh from `AbstractCellPopulation`
     * for use in this class
     */
    PottsMesh<DIM>* mpPottsMesh;

    /**
     * Pointer to a MutableMesh that can be created from the nodes of the PottsMesh in
     * order to solve PDEs on the population.
     */
    MutableMesh<DIM,DIM>* mpMutableMesh;

    /** The temperature of the system. Initialized to 0.1 in the constructor. */
    double mTemperature;

    /**
     * The number of MonteCarlo sweeps of the mesh performed each timestep.
     * Initialised to 1 in the constructor.
     */
    unsigned mNumSweepsPerTimestep;

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
        archive & boost::serialization::base_object<AbstractOnLatticeCellPopulation<DIM> >(*this);

        /*
         * In its current form the code does not allow the direct serialization
         * of the PottsMesh class, so instead we delete mpVoronoiTessellation.
         */
        delete mpElementTessellation;
        mpElementTessellation = nullptr;

        archive & mTemperature;
        archive & mNumSweepsPerTimestep;
    }

    /**
     * Check the consistency of internal data structures.
     * Each PottsElement must have a CellPtr associated with it.
     */
    void Validate();

    /**
     * Overridden WriteVtkResultsToFile() method.
     *
     * As a CellIdWriter is added to the population in the overridden method
     * OpenWritersFiles(), distinct cells can be visualized when viewing VTK
     * output of Potts model simulations.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     */
    virtual void WriteVtkResultsToFile(const std::string& rDirectory);

public:

    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely one CellPtr for each PottsElement in
     * the mesh.
     *
     * @param rMesh reference to a PottsMesh
     * @param rCells reference to a vector of CellPtrs
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     *                   (defaults to false)
     * @param validate whether to validate the cell population when it is created (defaults to true)
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    PottsBasedCellPopulation(PottsMesh<DIM>& rMesh,
                             std::vector<CellPtr>& rCells,
                             bool deleteMesh=false,
                             bool validate=true,
                             const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a vertex mesh.
     */
    PottsBasedCellPopulation(PottsMesh<DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~PottsBasedCellPopulation();

    /**
     * @return reference to mrMesh.
     */
    PottsMesh<DIM>& rGetMesh();

    /**
     * @return const reference to mrMesh (used in archiving).
     */
    const PottsMesh<DIM>& rGetMesh() const;

    /**
     * Overridden GetTetrahedralMeshForPdeModifier() method.
     *
     * @return a pointer to a tetrahedral mesh
     *
     * This method is called by AbstractGrowingDomainPdeModifier.
     */
    virtual TetrahedralMesh<DIM, DIM>* GetTetrahedralMeshForPdeModifier();

    /**
     * Get a particular PottsElement.
     *
     * @param elementIndex the global index of the PottsElement
     *
     * @return a pointer to the PottsElement.
     */
    PottsElement<DIM>* GetElement(unsigned elementIndex);

    /**
     * @return the number of PottsElements in the cell population.
     */
    unsigned GetNumElements();

    /**
     * Overridden GetNode() method.
     *
     * @param index global index of the specified node
     *
     * @return a pointer to the node.
     */
    Node<DIM>* GetNode(unsigned index);

    /**
     * Overridden GetNumNodes() method.
     *
     * @return the number of nodes in the cell population.
     */
    unsigned GetNumNodes();

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
     * Overridden GetLocationOfCellCentre() method.
     * Find where a given cell is in space.
     *
     * @param pCell the cell
     *
     * @return the location of the centre of mass of the element corresponding to this cell.
     */
    c_vector<double, DIM> GetLocationOfCellCentre(CellPtr pCell);

    /**
     * Get a pointer to the element corresponding to a given CellPtr.
     *
     * @param pCell the cell
     *
     * @return pointer to the element.
     */
    PottsElement<DIM>* GetElementCorrespondingToCell(CellPtr pCell);

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the cell population.
     *
     * @param pNewCell  the cell to add
     * @param pParentCell pointer to a parent cell (if required)
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    CellPtr AddCell(CellPtr pNewCell, CellPtr pParentCell=CellPtr());

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
     * Overridden UpdateCellLocations() method.
     *
     * @param dt time step
     */
    void UpdateCellLocations(double dt);

    /**
     * Overridden IsCellAssociatedWithADeletedLocation() method.
     *
     * @param pCell the cell
     * @return whether a given cell is associated with a deleted element.
     */
    bool IsCellAssociatedWithADeletedLocation(CellPtr pCell);

    /**
     * Remove the PottsElements which have been marked as deleted, and update the correspondence
     * with CellPtrs.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether cell population has had Births Or Deaths
     */
    void Update(bool hasHadBirthsOrDeaths=true);

    /**
     * Overridden OpenWritersFiles() method.
     *
     * Open all files in mCellPopulationWriters and mCellWriters for writing (not appending).
     *
     * @param rOutputFileHandler handler for the directory in which to open this file.
     */
    virtual void OpenWritersFiles(OutputFileHandler& rOutputFileHandler);

    /**
     * Overridden WriteResultsToFiles() method.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     */
    virtual void WriteResultsToFiles(const std::string& rDirectory);

    /**
     * A virtual method to accept a cell population writer so it can
     * write data from this object to file.
     *
     * @param pPopulationWriter the population writer.
     */
    virtual void AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter);

    /**
     * A virtual method to accept a cell population count writer so it can
     * write data from this object to file.
     *
     * @param pPopulationCountWriter the population count writer.
     */
    virtual void AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> > pPopulationCountWriter);

    /**
     * A virtual method to accept a cell writer so it can
     * write data from this object to file.
     *
     * @param pCellWriter the population writer.
     * @param pCell the cell whose data are being written.
     */
    virtual void AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell);

    /**
     * Overridden GetVolumeOfCell() method.
     *
     * @param pCell boost shared pointer to a cell
     * @return volume via associated mesh element
     */
    double GetVolumeOfCell(CellPtr pCell);

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
     * Outputs CellPopulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile);

    /**
     * Set mTemperature.
     *
     * @param temperature the temperature of the system
     */
    void SetTemperature(double temperature);

    /**
     * @return mTemperature
     */
    double GetTemperature();

    /**
     * Set mNumSweepsPerTimestep.
     *
     * @param numSweepsPerTimestep the number of MonteCarlo sweeps of the mesh performed each timestep
     */
    void SetNumSweepsPerTimestep(unsigned numSweepsPerTimestep);

    /**
     * @return mNumSweepsPerTimestep
     */
    unsigned GetNumSweepsPerTimestep();

    /**
     * Create a Element tessellation of the mesh for use in visualising the mesh.
     */
    void CreateElementTessellation();

    /**
     * @return a reference to mpElementTessellation.
     */
    VertexMesh<DIM,DIM>* GetElementTessellation();

    /**
     * Create a mutable mesh from the nodes of the Potts mesh. Used for solving PDEs
     */
    void CreateMutableMesh();

    /**
     * @return the Mutable Mesh
     */
    MutableMesh<DIM,DIM>* GetMutableMesh();

    /**
     * Overridden AddUpdateRule() method.
     *
     * @param pUpdateRule pointer to an update rule
     */
    virtual void AddUpdateRule(boost::shared_ptr<AbstractUpdateRule<DIM> > pUpdateRule);

    /**
     * Overridden GetCellDataItemAtPdeNode() method.
     *
     * @param pdeNodeIndex index of a node in a tetrahedral mesh for use
     *         with a PDE modifier
     * @param rVariableName the name of the cell data item to get
     * @param dirichletBoundaryConditionApplies where a Dirichlet boundary condition is used
     *        (optional; defaults to false)
     * @param dirichletBoundaryValue the value of the Dirichlet boundary condition, if used
     *        (optional; defaults to 0.0)
     *
     * @return the value of a CellData item (interpolated if necessary) at a node,
     *         specified by its index in a tetrahedral mesh for use with a PDE modifier.
     * This method can be called by PDE modifier classes.
     */
    virtual double GetCellDataItemAtPdeNode(unsigned pdeNodeIndex,
                                            std::string& rVariableName,
                                            bool dirichletBoundaryConditionApplies=false,
                                            double dirichletBoundaryValue=0.0);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PottsBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const PottsBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const PottsMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a PottsBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, PottsBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    PottsMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)PottsBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*POTTSBASEDCELLPOPULATION_HPP_*/
