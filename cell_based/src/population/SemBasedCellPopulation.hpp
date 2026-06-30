/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef SEMBASEDCELLPOPULATION_HPP_
#define SEMBASEDCELLPOPULATION_HPP_

#include "AbstractOffLatticeCellPopulation.hpp"
#include "SemMesh.hpp"
#include <string>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A facade class encapsulating a subcellular element method (SEM) based 
 * cell population.
 *
 * Contains a group of cells and maintains the associations between CellPtrs and 
 * SemElements in the SemMesh.
 */
template<unsigned DIM>
class SemBasedCellPopulation : public AbstractOffLatticeCellPopulation<DIM>
{
private:

    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this cell population has been de-serialized.
     */
    bool mDeleteMesh;

    /**
     * A static cast of the AbstractMesh from AbstractCellPopulation for use in 
     * this class.
     */
    SemMesh<DIM>* mpSemMesh;
    
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by 
     * load/save_construct_data.
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
    }

    /**
     * Check the consistency of internal cell-to-element maps.
     *
     * A SEM population owns one biological cell per live SemElement.
     * Validation checks that each live element has exactly one
     * associated CellPtr, that no cell maps to a deleted element, and that no
     * cell maps outside the mesh element range.
     */
    void Validate() override;

public:

    /**
     * Create a new cell population facade from a SemMesh and collection of 
     * CellPtrs.
     *
     * There must be precisely one CellPtr for each SemElement in the SemMesh.
     * If locationIndices is supplied, it defines the element index associated
     * with each cell; otherwise cells are attached to elements in order.
     *
     * @param rMesh reference to a SemMesh
     * @param rCells reference to a vector of CellPtrs
     * @param deleteMesh set to true if you want the cell population to free the 
     *                   mesh memory on destruction (defaults to false)
     * @param validate whether to validate the cell population when it is 
     *                 created (defaults to true)
     * @param locationIndices an optional vector of location indices that 
     *                        correspond to real cells
     */
    SemBasedCellPopulation(SemMesh<DIM>& rMesh,
                           std::vector<CellPtr>& rCells,
                           bool deleteMesh=false,
                           bool validate=true,
                           const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Constructor for use by boost serialization only.
     *
     * @param rMesh a SemMesh
     */
    SemBasedCellPopulation(SemMesh<DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~SemBasedCellPopulation();

    /**
     * @return reference to the underlying SEM mesh.
     */
    SemMesh<DIM>& rGetMesh();

    /**
     * @return const reference to the underlying SEM mesh, used when archiving.
     */
    const SemMesh<DIM>& rGetMesh() const;

    /**
     * Get a particular SemElement.
     *
     * @param elementIndex the global index of the SemElement
     *
     * @return a pointer to the SemElement.
     */
    SemElement<DIM>* GetElement(unsigned elementIndex);

    /**
     * Overridden GetNumNodes() method.
     *
     * @return the number of nodes in the cell population.
     */
    unsigned GetNumNodes() override;
    
    /**
     * GetNumElements() method.
     *
     * @return the number of SemElements in the cell population.
     */
    unsigned GetNumElements();

    /**
     * Overridden GetLocationOfCellCentre() method.
     *
     * Find the centre of mass of a given cell (assuming uniform density).
     * Note that, as there is no guarantee of convexity, this may lie
     * outside the SemElement corresponding to the cell.
     *
     * @param pCell pointer to a cell in the population
     *
     * @return the location of the centre of mass of the SemElement 
     *         corresponding to this cell.
     */
    c_vector<double, DIM> GetLocationOfCellCentre(CellPtr pCell) override;

    /**
     * Overridden GetNode() method.
     *
     * @param index global index of the specified node
     *
     * @return a pointer to the node.
     */
    Node<DIM>* GetNode(unsigned index) override;

    /**
     * Overridden GetNeighbouringLocationIndices() method.
     *
     * Given a cell, returns the set of element location indices corresponding
     * to neighbouring cells. SEM elements are treated as neighbours when
     * their nodes appear in the current node-pair interaction
     * list, so Update() should be called before this method is used for force
     * or topology queries.
     *
     * @param pCell pointer to a cell in the population
     * @return the set of neighbouring location indices.
     */
    std::set<unsigned> GetNeighbouringLocationIndices(CellPtr pCell) override;

    /**
     * Overridden AddNode() method.
     *
     * Add a new node to the cell population.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in cell population
     */
    unsigned AddNode(Node<DIM>* pNewNode) override;

    /**
     * Overridden SetNode() method.
     *
     * Move the node with a given index to a new point in space.
     *
     * @param index the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned index, ChastePoint<DIM>& rNewLocation) override;

    /**
     * Get a pointer to the SemElement corresponding to a given CellPtr.
     *
     * @param pCell pointer to a cell in the population
     *
     * @return pointer to the SemElement.
     */
    SemElement<DIM>* GetElementCorrespondingToCell(CellPtr pCell);

    /**
     * Overridden GetVolumeOfCell() method.
     *
     * @param pCell pointer to a cell in the population
     * 
     * @return volume via associated SemElement.
     */
    double GetVolumeOfCell(CellPtr pCell) override;

    /**
     * Overridden OutputCellPopulationParameters() method.
     *
     * Writes SEM population parameters and then delegates shared off-lattice
     * parameters, such as damping constants, to the base class.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile) override;

    /**
     * Overridden WriteVtkResultsToFile() method.
     *
     * Writes the current SEM mesh as VTK point/element data for a
     * simulation time step, including any registered node point data writers,
     * and appends the generated file to the VTK meta-file.
     *
     * @param rDirectory output directory for the VTK files
     */
    void WriteVtkResultsToFile(const std::string& rDirectory) override;

    /**
     * Overridden GetTetrahedralMeshForPdeModifier() method.
     *
     * SEM populations do not currently expose a tetrahedral mesh for PDE
     * coupling, so this returns nullptr.
     *
     * @return a pointer to a PDE mesh, or nullptr when unsupported
     */
    TetrahedralMesh<DIM, DIM>* GetTetrahedralMeshForPdeModifier() override;

    /**
     * Overridden GetCellDataItemAtPdeNode() method.
     *
     * SEM PDE coupling is not currently implemented, so this returns a neutral
     * value rather than interpolating cell data onto PDE nodes.
     *
     * @param pdeNodeIndex index of the PDE node
     * @param item cell data item name
     * @param averageOverNeighbours whether to average over neighbouring cells
     * @param defaultValue default value for missing data
     * @return the requested cell data value at the PDE node
     */
    double GetCellDataItemAtPdeNode(unsigned pdeNodeIndex, std::string& item, bool averageOverNeighbours, double defaultValue) override;

    /**
     * Overridden IsCellAssociatedWithADeletedLocation() method.
     *
     * Checks whether the SemElement representing a cell has been
     * marked deleted during population cleanup.
     *
     * @param pCell pointer to a cell in the population
     * @return whether the cell maps to a deleted SemElement
     */
    bool IsCellAssociatedWithADeletedLocation(CellPtr pCell) override;

    /**
     * Overridden AddCell() method.
     *
     * This would add a new biological cell by creating or
     * dividing an SemElement. SEM element division is not currently implemented,
     * so the method throws rather than creating an inconsistent map-only cell.
     *
     * @param pNewCell pointer to the new cell
     * @param pParentCell pointer to the parent cell
     * @return pointer to the created cell if addition is supported
     */
    CellPtr AddCell(CellPtr pNewCell, CellPtr pParentCell) override;

    /**
     * Overridden GetDefaultTimeStep() method.
     *
     * @return the default mechanics time step for SEM simulations.
     */
    double GetDefaultTimeStep() override;

    /**
     * Overridden RemoveDeadCells() method.
     *
     * Removes biological cells that have died, marks their
     * associated SemElements deleted, unregisters those elements from their
     * nodes, and removes the cell-to-element map entries.
     *
     * @return the number of cells removed
     */
    unsigned RemoveDeadCells() override;

    /**
     * Overridden Update() method.
     *
     * Refreshes spatial search data for the current node
     * positions, then rebuilds the node-pair interaction list used by SEM
     * force and neighbour queries.
     *
     * @param hasHadBirthsOrDeaths whether the population changed size this step
     */
    void Update(bool hasHadBirthsOrDeaths) override;

    /**
     * Overridden GetWidth() method.
     *
     * Delegates to the underlying mesh to measure the population extent in a
     * coordinate direction.
     *
     * @param rDimension dimension over which to measure width
     * @return width of the SEM mesh in that dimension
     */
    double GetWidth(const unsigned& rDimension) override;

    /**
     * Overridden GetNeighbouringNodeIndices() method.
     *
     * Returns the nodes that are currently close enough to
     * interact with the supplied node, using the node-pair list rebuilt by
     * Update().
     *
     * @param index global node index
     * @return set of neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index) override;

    /**
     * Accept a population writer.
     *
     * This is the visitor dispatch point used by output writers
     * that operate on the whole SEM population.
     *
     * @param pPopulationWriter population writer to visit this population
     */
    void AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter) override;

    /**
     * Accept a population count writer.
     *
     * This is the visitor dispatch point used by writers that
     * count cells or cell properties across the SEM population.
     *
     * @param pPopulationCountWriter population count writer to visit this population
     */
    void AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> > pPopulationCountWriter) override;

    /**
     * Accept a population event writer.
     *
     * This is the visitor dispatch point used by writers that
     * record population-level cell events.
     *
     * @param pPopulationEventWriter population event writer to visit this population
     */
    void AcceptPopulationEventWriter(boost::shared_ptr<AbstractCellPopulationEventWriter<DIM, DIM> > pPopulationEventWriter) override;

    /**
     * Accept a cell writer for a single cell.
     *
     * This dispatches cell-level output for a cell while providing
     * the SEM population context needed to interpret its location.
     *
     * @param pCellWriter cell writer to visit the cell
     * @param pCell cell to write
     */
    void AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell) override;

    /**
     * Overridden CheckForStepSizeException() method.
     *
     * This is the hook where SEM-specific movement restrictions
     * would reject a node displacement that is too large for the current time
     * step. No SEM-specific restriction is currently applied here.
     *
     * @param nodeIndex index of the node being moved
     * @param rDisplacement proposed node displacement
     * @param dt current time step
     */
    void CheckForStepSizeException(unsigned nodeIndex, c_vector<double,DIM>& rDisplacement, double dt) override;

    /**
     * Overridden GetDampingConstant() method.
     *
     * Returns the drag coefficient used when updating a SEM node.
     * For nodes contained in multiple live elements, the damping constant is
     * averaged over the cells associated with those elements, using mutant or
     * normal damping according to each cell mutation state.
     *
     * @param nodeIndex global node index
     * @return damping constant for the node
     */
    double GetDampingConstant(unsigned nodeIndex) override;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SemBasedCellPopulation.
 *
 * The population facade is reconstructed around its mesh on load, so the
 * archive stores the mesh pointer needed by load_construct_data().
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SemBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const SemMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a SemBasedCellPopulation.
 *
 * This retrieves the archived SEM mesh and invokes the serialization-only
 * constructor in-place, after which Boost restores the base class cell and
 * mapping state.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SemBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    SemMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)SemBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*SEMBASEDCELLPOPULATION_HPP_*/
