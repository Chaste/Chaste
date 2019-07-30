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

#ifndef CABASEDCELLPOPULATION_HPP_
#define CABASEDCELLPOPULATION_HPP_

#include "AbstractOnLatticeCellPopulation.hpp"
#include "PottsMesh.hpp"
#include "VertexMesh.hpp"
#include "AbstractUpdateRule.hpp"
#include "AbstractCaBasedDivisionRule.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template<unsigned DIM> class AbstractCaBasedDivisionRule; // Circular definition thing.

/**
 * A facade class encapsulating a cell population under the Cellular
 * Automaton (CA) framework.
 *
 * Contains a group of cells and maintains the associations
 * between CellPtrs and nodes in a specialised PottsMesh class.
 *
 * When used here the PottsMesh has no elements as Cells are associated with nodes.
 * The PottsMesh is used to define node connectivity.
 *
 * Multiple cells can be associated at a single node.
 */
template<unsigned DIM>
class CaBasedCellPopulation : public AbstractOnLatticeCellPopulation<DIM>
{
    friend class TestCaBasedCellPopulation;

private:

    /** The carrying capacity (number of cells allowed per site). */
    unsigned mLatticeCarryingCapacity;

    /**
     * The update rules used to determine the new location of the cells.
     * These rules specify is cells switch locations.
     */
    std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > > mSwitchingUpdateRuleCollection;

    /** Records for each node the node the number of spaces available. */
    std::vector<unsigned> mAvailableSpaces;

    /** A pointer to a division rule that is used to specify how cells divide. I.e do they move other cells out of the way.
     * This is a specialisation for CA models. */
    boost::shared_ptr<AbstractCaBasedDivisionRule<DIM> > mpCaBasedDivisionRule;

    /**
     * Set the empty sites by taking in a set of which nodes indices are empty sites.
     *
     * @param rEmptySiteIndices set of node indices corresponding to empty sites
     */
    void SetEmptySites(const std::set<unsigned>& rEmptySiteIndices);

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
        archive & mSwitchingUpdateRuleCollection;
        archive & mLatticeCarryingCapacity;
        archive & mAvailableSpaces;
        archive & mpCaBasedDivisionRule;
    }

    /**
     * Overridden Validate() method.
     *
     * Not used in CA simulations so just contains NEVER_REACHED
     */
    void Validate();

    /**
     * Overridden WriteVtkResultsToFile() method.
     *
     * This method offsets cells so can visulaise multiple cells at a single site
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     */
    virtual void WriteVtkResultsToFile(const std::string& rDirectory);

public:

    /**
     * Create a new cell population facade from a mesh, a vector of location indices
     * and a collection of cells.
     *
     * There must be precisely one CellPtr for each entry of the locationIndices vector.
     *
     * @param rMesh reference to a PottsMesh
     * @param rCells reference to a vector of CellPtrs
     * @param locationIndices a vector of location indices that correspond to real cells
     * @param latticeCarryingCapacity an optional parameter to allow more than one cell per site
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     *                   (defaults to false)
     * @param validate whether to validate the cell population when it is created (defaults to false as not used in CA simulations)
     */
    CaBasedCellPopulation(PottsMesh<DIM>& rMesh,
                                  std::vector<CellPtr>& rCells,
                                  const std::vector<unsigned> locationIndices,
                                  unsigned latticeCarryingCapacity=1u,
                                  bool deleteMesh=false,
                                  bool validate=false);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a Ca mesh.
     */
    CaBasedCellPopulation(PottsMesh<DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~CaBasedCellPopulation();

    /**
     * @return mAvailableSpaces.
     */
    std::vector<unsigned>& rGetAvailableSpaces();

    /**
     * Find if a given node has space available. Overridden in subclasses to implement different division regimes.
     *
     * \todo Which subclasses? Why is the second input argument needed?
     *
     * @param index  The global index of a specified node.
     * @param pCell  The cell wanting to divide into the lattice site (defaults to NULL).
     *
     * @return whether the node is an empty site
     */
    virtual bool IsSiteAvailable(unsigned index, CellPtr pCell);

    /**
     * @return reference to #mrMesh.
     */
    PottsMesh<DIM>& rGetMesh();

    /**
     * @return const reference to #mrMesh (used in archiving).
     */
    const PottsMesh<DIM>& rGetMesh() const;

    /**
     * Overridden GetTetrahedralMeshForPdeModifier() method.
     *
     * @return a pointer to a tetrahedral mesh, for use with a PDE modifier.
     *
     * This method is called by AbstractGrowingDomainPdeModifier.
     */
    virtual TetrahedralMesh<DIM, DIM>* GetTetrahedralMeshForPdeModifier();

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
     * Note: In keeping with other parts of the code in this class, we assume a Moore
     * neighbourhood. Also, at present this method assumes a unit carrying capacity at
     * each lattice site.
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
     * @return the location of the node corresponding to this cell.
     */
    c_vector<double, DIM> GetLocationOfCellCentre(CellPtr pCell);

    /**
     * Overridden AddCellUsingLocationIndex method to add a cell to a given location index.
     * Also updates mAvailableSpaces
     *
     * @param index the location index
     * @param pCell the cell.
     */
    void AddCellUsingLocationIndex(unsigned index, CellPtr pCell);

    /**
     * Overridden AddCellUsingLocationIndex method to remove a cell from a given location index.
     * Also updates mAvailableSpaces
     *
     * @param index the location index
     * @param pCell the cell.
     */
    void RemoveCellUsingLocationIndex(unsigned index, CellPtr pCell);

    /**
     * Get a pointer to the node corresponding to a given CellPtr.
     *
     * @param pCell the cell
     *
     * @return pointer to the node.
     */
     Node<DIM>* GetNodeCorrespondingToCell(CellPtr pCell);

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
      * Calculate the propensity of a dividing into a given site.
      * Overridden in child classes to define other division methods, e.g. directed division.
      *
      * \todo This functionality should be moved into the CA-based division rule hierarchy
      *
      * @param currentNodeIndex The index of the current node/lattice site
      * @param targetNodeIndex The index of the target node/lattice site
      * @param pCell a pointer to the cell (needed if more than one cell per lattice site
      * @return The probability of the cell dividing from the current node to the target node
      */
     double virtual EvaluateDivisionPropensity(unsigned currentNodeIndex,
                                               unsigned targetNodeIndex,
                                               CellPtr pCell);
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
     * Overridden OpenWritersFiles() method.
     *
     * Open all files in mCellPopulationWriters and mCellWriters for writing (not appending).
     *
     * @param rOutputFileHandler handler for the directory in which to open this file.
     */
    virtual void OpenWritersFiles(OutputFileHandler& rOutputFileHandler);

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
     * @return whether a given cell is associated with a deleted node.
     */
    bool IsCellAssociatedWithADeletedLocation(CellPtr pCell);

    /**
     * Overridden Update() method.
     *
     * Checks association of nodes with CellPtrs.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether cell population has had Births Or Deaths
     */
    void Update(bool hasHadBirthsOrDeaths=true);

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
     * Note this returns the size of the underlying mesh not the population of cells so here it will be the same for all time.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    double GetWidth(const unsigned& rDimension);

    /**
     * Overridden RemoveAllUpdateRules() method.
     *
     * Remove any update rules previously passed to this population by clearing
     * mUpdateRuleCollection and mSwitchingUpdateRuleCollection.
     */
    void RemoveAllUpdateRules();

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
     * Overridden IsRoomToDivide() method.
     * Returns whether there are any available neighbouring sites to the
     * one occupied by a given cell.
     *
     * @param pCell pointer to a cell
     * @return whether the cell has any free neighbouring sites
     */
    bool IsRoomToDivide(CellPtr pCell);

    /**
     * @return The Ca division rule that is currently being used.
     */
    boost::shared_ptr<AbstractCaBasedDivisionRule<DIM> > GetCaBasedDivisionRule();

    /**
     * Set the division rule for this population.
     *
     * @param pCaBasedDivisionRule  pointer to the new division rule
     */
    void SetCaBasedDivisionRule(boost::shared_ptr<AbstractCaBasedDivisionRule<DIM> > pCaBasedDivisionRule);

    /**
     * Overridden AddUpdateRule() method.
     *
     * @param pUpdateRule pointer to an update rule
     */
    virtual void AddUpdateRule(boost::shared_ptr<AbstractUpdateRule<DIM> > pUpdateRule);

    /**
     * Overridden AddUpdateRule() method.
     *
     * Get the collection of update rules to be used with this population.
     * This vector is comprised of mUpdateRuleCollection and mSwitchingUpdateRuleCollection,
     * one after the other.
     *
     * @return the update rule collection
     */
    virtual const std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > > GetUpdateRuleCollection() const;

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

    /**
     * Overridden IsPdeNodeAssociatedWithNonApoptoticCell() method.
     *
     * @param pdeNodeIndex index of a node in a tetrahedral mesh for use with a PDE modifier
     *
     * @return if a node, specified by its index in a tetrahedral mesh for use
     *         with a PDE modifier, is associated with a non-apoptotic cell.
     * This method can be called by PDE classes.
     */
    virtual bool IsPdeNodeAssociatedWithNonApoptoticCell(unsigned pdeNodeIndex);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CaBasedCellPopulation)


namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CaBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CaBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const PottsMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a CaBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CaBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    PottsMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)CaBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...


#endif /*CABASEDCELLPOPULATION_HPP_*/
