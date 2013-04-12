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

#ifndef ABSTRACTCELLPOPULATION_HPP_
#define ABSTRACTCELLPOPULATION_HPP_

#include "Cell.hpp"
#include "OutputFileHandler.hpp"

#include <list>
#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "Node.hpp"
#include "CellPropertyRegistry.hpp"
#include "Identifiable.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ApoptoticCellProperty.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"
#include "AbstractMesh.hpp"

/** Forward declaration of cell population writers. */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCellPopulationWriter;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCellWriter;
/**
 * An abstract facade class encapsulating a cell population.
 *
 * Contains a group of cells and associated methods.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractCellPopulation : public Identifiable
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCells;
        archive & mLocationCellMap;
        archive & mCellLocationMap;
        archive & mCellCyclePhaseCount;
        archive & mpCellPropertyRegistry;
        archive & mOutputResultsForChasteVisualizer;
        archive & mOutputCellIdData;
        archive & mOutputCellMutationStates;
        archive & mOutputCellAncestors;
        archive & mOutputCellProliferativeTypes;
        archive & mOutputCellVariables;
        archive & mOutputCellCyclePhases;
        archive & mOutputCellAges;
        archive & mOutputCellVolumes;
    }

protected:

    /** Map location (node or VertexElement) indices back to cells */
    std::map<unsigned, std::set<CellPtr> > mLocationCellMap;

    /** Map cells to location (node or VertexElement) indices */
    std::map<Cell*, unsigned> mCellLocationMap;

    /** Reference to the mesh */
    AbstractMesh<ELEMENT_DIM, SPACE_DIM>& mrMesh;

    /** List of cells */
    std::list<CellPtr> mCells;

    /** Current cell cycle phase counts */
    std::vector<unsigned> mCellCyclePhaseCount;

    /** Population centroid */
    c_vector<double, SPACE_DIM> mCentroid;

    /** Results file for node visualization */
    out_stream mpVizNodesFile;

    /** Results file for cell visualization */
    out_stream mpVizCellProliferativeTypesFile;

    /** Results file for cell visualization */
    out_stream mpVizCellProliferativePhasesFile;

    /** Results file for cell mutation states */
    out_stream mpCellMutationStatesFile;

    /** Results file for cell ancestors */
    out_stream mpVizCellAncestorsFile;

    /** Results file for cell types */
    out_stream mpCellProliferativeTypesFile;

    /** Results file for cell cycle phases */
    out_stream mpCellCyclePhasesFile;

    /** Results file for cell variables */
    out_stream mpCellVariablesFile;

    /** Results file for cell ages */
    out_stream mpCellAgesFile;

    /** Results file for logged cell data. */
    out_stream mpCellIdFile;

    /** Results file for cell volume (in 3D) or area (in 2D) data. */
    out_stream mpCellVolumesFile;

    /** Results file for boundary nodes. */
    out_stream mpVizBoundaryNodesFile;

    /** A cache of where the results are going (used for VTK writer). */
    std::string mDirPath;

    /** Meta results file for VTK. */
    out_stream mpVtkMetaFile;

    /** Cell property registry. */
    boost::shared_ptr<CellPropertyRegistry> mpCellPropertyRegistry;

    /** Whether to write results to file for visualization using the Chaste java visualizer (defaults to true). */
    bool mOutputResultsForChasteVisualizer;

    /** Whether to write cell ID data to file. */
    bool mOutputCellIdData;

    /** Whether to count the number of each cell mutation state and output to file. */
    bool mOutputCellMutationStates;

    /** Whether to output the ancestor of each cell to a visualizer file. */
    bool mOutputCellAncestors;

    /** Whether to count the number of each cell type and output to file. */
    bool mOutputCellProliferativeTypes;

    /** Whether to write the cell variables to a file. */
    bool mOutputCellVariables;

    /** Whether to write the cell cycle phases to a file. */
    bool mOutputCellCyclePhases;

    /** Whether to write the cell ages to a file. */
    bool mOutputCellAges;

    /** Whether to write the cell volumes (in 3D) or areas (in 2D) to a file. */
    bool mOutputCellVolumes;

    /**
     * Check consistency of our internal data structures.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void Validate()=0;

    /**
     * Write the current results to mpVtkMetaFile.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void WriteVtkResultsToFile()=0;

    /**
     * Constructor that just takes in a mesh.
     *
     * @param rMesh the mesh for the population.
     */
    AbstractCellPopulation(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

public:

    /**
     * AbstractCellPopulation Constructor.
     *
     * @note Warning: the passed-in vector of cells will be emptied, even if the constructor
     * throws an exception!
     *
     * @param rMesh a refernce to the mesh underlying the cell population
     * @param rCells a vector of cells.  Copies of the cells will be stored in the cell population,
     *     and the passed-in vector cleared.
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    AbstractCellPopulation( AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                            std::vector<CellPtr>& rCells,
                            const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractCellPopulation();

    /**
     * Initialise each cell's cell-cycle model.
     */
    void InitialiseCells();

    /**
     * Add an item of cell data to every cell in the population
     * @param dataName is the name associated with the data
     * @param dataValue is the value of the data, initially the same for each cell
     */
    void SetDataOnAllCells(const std::string& dataName, double dataValue);

    /**
     * @return reference to the mesh, mrMesh.
     */
    AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rGetMesh();

    /**
     * @return reference to mCells.
     */
    std::list<CellPtr>& rGetCells();

    /**
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @return the number of nodes in the cell population.
     */
    virtual unsigned GetNumNodes()=0;

    /**
     * Find where a given cell is in space.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCell the cell
     * @return the location of the cell
     */
    virtual c_vector<double, SPACE_DIM> GetLocationOfCellCentre(CellPtr pCell)=0;

    /**
     * Get a pointer to the node with a given index.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param index  global index of the specified node
     * @return a pointer to the node with a given index.
     */
    virtual Node<SPACE_DIM>* GetNode(unsigned index)=0;

    /**
     * Move the node with a given index to a new point in space.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    virtual void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM>& rNewLocation)=0;

    /**
     * Helper method for establishing if a cell is real.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCell the cell
     * @return whether a given cell is associated with a deleted
     *         node (cell-centre models) or element (vertex models).
     */
    virtual bool IsCellAssociatedWithADeletedLocation(CellPtr pCell)=0;

    /**
     * Add a new cell to the cell population.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pNewCell  the cell to add
     * @param rCellDivisionVector  a vector providing information regarding how the cell division should occur
     *     (for cell-centre cell populations, this vector is the position of the daughter cell; for vertex cell populations it
     *      can be used by any subclass of CellBasedSimulation to as a means of dictating the axis along which
     *      the parent cell divides)
     * @param pParentCell pointer to a parent cell (if required)
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way).
     */
    virtual CellPtr AddCell(CellPtr pNewCell, const c_vector<double,SPACE_DIM>& rCellDivisionVector, CellPtr pParentCell=CellPtr())=0;

    class Iterator; // Forward declaration; see below

    /**
     * Remove all cells labelled as dead.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @return number of cells removed
     */
    virtual unsigned RemoveDeadCells()=0;

    /**
     * Remove the Nodes (for cell-centre) or VertexElements (for cell-vertex) which
     * have been marked as deleted and update the correspondence with Cells.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether cell population has had Births Or Deaths
     */
    virtual void Update(bool hasHadBirthsOrDeaths=true)=0;

    /**
     * Find out how many cells of each mutation state there are
     *
     * @return The number of cells of each mutation state (evaluated at each visualizer output), with default ordering
     * [0] = healthy count
     * [1] = APC one hit
     * [2] = APC two hit
     * [3] = beta catenin one hit
     */
    std::vector<unsigned> GetCellMutationStateCount();

    /**
     * Find out how many cells of each type there are.
     *
     * @return The number of cells of each type (evaluated at each visualizer output), with default ordering
     * [0] = STEM
     * [1] = TRANSIT
     * [2] = DIFFERENTIATED
     * [3] = DEFAULT
     */
    std::vector<unsigned> GetCellProliferativeTypeCount();

    /**
     * Find out how many cells in each cell cycle phase there are.
     *
     * @return The number of cells of each phase (evaluated at each visualizer output)
     * [0] = G_ZERO_PHASE
     * [1] = G_ONE_PHASE
     * [2] = S_PHASE
     * [3] = G_TWO_PHASE
     * [4] = M_PHASE
     */
    const std::vector<unsigned>& rGetCellCyclePhaseCount() const;

    /**
     * @return the number of real cells.
     */
    unsigned GetNumRealCells();

    /**
     * Sets the Ancestor index of all the cells at this time to be the
     * same as their location index, can be used to trace clonal populations.
     */
    void SetCellAncestorsToLocationIndices();

    /**
     * Write cell ID data to mpCellIdFile.
     */
    void WriteCellIdDataToFile();

    /**
     * Loops over cells and makes a list of the ancestors that
     * are part of the cell population.
     *
     * @return remaining_ancestors  The size of this set tells you how many clonal populations remain.
     */
    std::set<unsigned> GetCellAncestors();

    /**
     * Get the cell corresponding to a given location index.
     *
     * This method assumes that there is at most one cell attached to a location index and an assertion fails if not.
     * \todo should be an exception?
     *
     * @param index the location index
     *
     * @return the cell.
     */
    CellPtr GetCellUsingLocationIndex(unsigned index);

    /**
     * Get the set of cells corresponding to a given location index.
     *
     * Note the set may be empty.
     *
     * @param index the location index
     *
     * @return the set of cells.
     */
    std::set<CellPtr> GetCellsUsingLocationIndex(unsigned index);

    /**
     * Returns whether or not a cell is associated with a location index
     *
     * @param index the location index
     *
     * @return the cell.
     */
    bool IsCellAttachedToLocationIndex(unsigned index);

    /**
     * Set the cell corresponding to a given location index.
     *
     * Assumes there is one cell for each location index and replaces any existing cell attached to the location index.
     * If you wish to attach an additional cell to a location index use AddCellUsingLocaitonIndex as SetCellUsingLocation
     * Index will overwrite cells attached to this index.
     *
     * @param index the location index
     * @param pCell the cell.
     */
    void SetCellUsingLocationIndex(unsigned index, CellPtr pCell);

    /**
     * Adds a cell to a given location index.
     *
     * @param index the location index
     * @param pCell the cell.
     */
    virtual void AddCellUsingLocationIndex(unsigned index, CellPtr pCell);

    /**
     * Removes a cell from a given location index.
     *
     * @param index the location index
     * @param pCell the cell.
     */
    virtual void RemoveCellUsingLocationIndex(unsigned index, CellPtr pCell);

    /**
     * Change the location index of a cell in mLocationCellMap and mCellLocationMap
     *
     * @param pCell the cell to move
     * @param old_index the old location index
     * @param new_index the new location index
     */
    void MoveCellInLocationMap(CellPtr pCell, unsigned old_index, unsigned new_index);

    /**
     * Get the location index corresponding to a given cell.
     *
     * Assumes there is one location index for each cell and an assertion fails if not.
     *
     * @param pCell the cell
     *
     * @return the location index.
     */
    unsigned GetLocationIndexUsingCell(CellPtr pCell);

    /**
     * @return registry of cell properties used in this cell population.
     */
    boost::shared_ptr<CellPropertyRegistry> GetCellPropertyRegistry();

    /**
     * Set a default ordering on cell mutation states and cell proliferative types, so that
     * existing tests don't need to specify the old ordering explicitly.
     */
    void SetDefaultCellMutationStateAndProliferativeTypeOrdering();

    /**
     * Calculate the 'width' of any dimension of the cell population.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    virtual double GetWidth(const unsigned& rDimension)=0;

    /**
     * @return the volume (or area in 2D, or length in 1D) of a given cell.
     *
     * As this method is pure virtual, it must be overridden in subclasses.
     *
     * @param pCell boost shared pointer to a cell
     */
    virtual double GetVolumeOfCell(CellPtr pCell)=0;

    /**
     * Given a node index, returns the set of neighbouring node indices.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param index the node index
     * @return the set of neighbouring node indices.
     */
    virtual std::set<unsigned> GetNeighbouringNodeIndices(unsigned index)=0;

    /**
     * @return the centroid of the cell population.
     */
    c_vector<double, SPACE_DIM> GetCentroidOfCellPopulation();

    /**
     * Update the ownership of cell in a parallel cell-based simulation.
     */
    virtual void UpdateCellProcessLocation();

    /**
     * Use an output file handler to create output files for visualizer and post-processing.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     * @param cleanOutputDirectory  whether to delete the contents of the output directory prior to output file creation
     */
    virtual void CreateOutputFiles(const std::string& rDirectory,
                                   bool cleanOutputDirectory);

    /**
     * Clear the counters used for cell population output.
     */
    void ResetCellCounters();

    /**
     * Write results from the current cell population state to output files.
     */
    virtual void WriteResultsToFiles();

    /**
     * A virtual method to accept a cell population writer so it can
     * write data from this object to file.
     *
     * @param pPopulationWriter the population writer.
     */
    virtual void AcceptPopulationWriter(AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>* pPopulationWriter)=0;

    /**
     * A virtual method to accept a cell writer so it can
     * write data from this object to file.
     *
     * @param pCellWriter the population writer.
     */
    virtual void AcceptCellWriter(AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>* pCellWriter)=0;

    /**
     * Write the current time and node results to output files.
     */
    virtual void WriteTimeAndNodeResultsToFiles();

    /**
     * Calls GenerateCellResults() on each cell then calls WriteCellResultsToFiles().
     */
    virtual void GenerateCellResultsAndWriteToFiles()=0;

    /**
     * Generate results for a given cell in the current cell population state to output files.
     *
     * @param pCell pointer to the cell
     */
    virtual void GenerateCellResults(CellPtr pCell);

    /**
     * Write the current volume of each cell to file.
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void WriteCellVolumeResultsToFile()=0;

    /**
     * Write the current state of each cell to output files.
     *
     */
    void WriteCellResultsToFiles();

    /**
     * Close any output files.
     */
    virtual void CloseOutputFiles();

    /**
     * Outputs CellPopulation used in the simulation to file and then calls OutputCellPopulationParameters to output all relevant parameters.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationInfo(out_stream& rParamsFile);

    /**
     * Outputs CellPopulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationParameters(out_stream& rParamsFile)=0;

    /**
     * @return mOutputResultsForChasteVisualizer
     */
    bool GetOutputResultsForChasteVisualizer();

    /**
     * @return mOutputCellIdData
     */
    bool GetOutputCellIdData();

    /**
     * @return mOutputCellMutationStates
     */
    bool GetOutputCellMutationStates();

    /**
     * @return mOutputCellAncestors
     */
    bool GetOutputCellAncestors();

    /**
     * @return mOutputCellProliferativeTypes
     */
    bool GetOutputCellProliferativeTypes();

    /**
     * @return mOutputCellVariables
     */
    bool GetOutputCellVariables();

    /**
     * @return mOutputCellCyclePhases
     */
    bool GetOutputCellCyclePhases();

    /**
     * @return mOutputCellAges
     */
    bool GetOutputCellAges();

    /**
     * @return mOutputCellVolumes
     */
    bool GetOutputCellVolumes();

    /**
     * Set mOutputResultsForChasteVisualizer.
     *
     * @param outputResultsForChasteVisualizer the new value of mOutputResultsForChasteVisualizer
     */
    void SetOutputResultsForChasteVisualizer(bool outputResultsForChasteVisualizer);

    /**
     * Set mOutputCellIdData.
     *
     * @param outputCellIdData the new value of mOutputCellIdData
     */
    void SetOutputCellIdData(bool outputCellIdData);

    /**
     * Set mOutputCellMutationStates.
     *
     * @param outputCellMutationStates the new value of mOutputCellMutationStates
     */
    void SetOutputCellMutationStates(bool outputCellMutationStates);

    /**
     * Set mOutputCellAncestors.
     *
     * @param outputCellAncestors the new value of mOutputCellAncestors
     */
    void SetOutputCellAncestors(bool outputCellAncestors);

    /**
     * Set mOutputCellProliferativeTypes.
     *
     * @param outputCellProliferativeTypes the new value of mOutputCellProliferativeTypes
     */
    void SetOutputCellProliferativeTypes(bool outputCellProliferativeTypes);

    /**
     * Set mOutputCellVariables.
     *
     * @param outputCellVariables the new value of mOutputCellVariables
     */
    void SetOutputCellVariables(bool outputCellVariables);

    /**
     * Set mOutputCellCyclePhases.
     *
     * @param outputCellCyclePhases the new value of mOutputCellCyclePhases
     */
    void SetOutputCellCyclePhases(bool outputCellCyclePhases);

    /**
     * Set mOutputCellAges.
     *
     * @param outputCellAges the new value of mOutputCellAges
     */
    void SetOutputCellAges(bool outputCellAges);

    /**
     * Set mOutputCellVolumes.
     *
     * @param outputCellVolumes the new value of mOutputCellVolumes
     */
    void SetOutputCellVolumes(bool outputCellVolumes);

    /**
     * @return The width (maximum distance to centroid) of the cell population
     *     in each dimension
     */
    c_vector<double,SPACE_DIM> GetSizeOfCellPopulation();

    /**
     * @return whether there is room into which a given cell may divide.
     * Returns true by default, but may be overridden in subclasses.
     *
     * @param pCell pointer to a cell
     */
    virtual bool IsRoomToDivide(CellPtr pCell);

    /**
     * Iterator class allows one to iterate over cells in the cell population.
     * Dereferencing the iterator will give you the current cell.
     * There are also methods to get the node representing this cell,
     * and the location of that node.
     */
    class Iterator
    {
    public:

        /**
         * Dereference the iterator giving you a pointer to the current cell.
         * @return reference
         */
        inline CellPtr operator*();

        /**
         * Unusually for an iterator over a collection of pointers, this method
         * allows you to access the object pointed at, rather than the pointer itself.
         * @return pointer
         */
        inline CellPtr operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         * @return not-equal
         */
        inline bool operator!=(const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator& rOther);

        /**
         * Prefix increment operator.
         * @return incremented object
         */
        inline Iterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * @param rCellPopulation the cell population
         * @param cellIter iterator over list of cells
         */
        Iterator(AbstractCellPopulation& rCellPopulation, std::list<CellPtr>::iterator cellIter);

        /**
         * The iterator must have a virtual destructor.
         */
        virtual ~Iterator()
        {}

    private:

        /**
         * Private helper function.
         * @return whether we're pointing at a real cell.
         * Assumes we are within range (i.e. not at End).
         *
         * Real cells are not deleted.
         */
        virtual inline bool IsRealCell();

        /**
         * Private helper function.
         * @return whether we're at the end of the cells.
         */
        inline bool IsAtEnd();

        /** The cell population member. */
        AbstractCellPopulation& mrCellPopulation;

        /** Cell iterator member. */
        std::list<CellPtr>::iterator mCellIter;
    };

    /**
     * @return iterator pointing to the first cell in the cell population
     */
    Iterator Begin();

    /**
     * @return iterator pointing to one past the last cell in the cell population
     */
    Iterator End();
};

enum cell_colours
{
    STEM_COLOUR, // 0
    TRANSIT_COLOUR, // 1
    DIFFERENTIATED_COLOUR, // 2
    EARLY_CANCER_COLOUR, // 3
    LATE_CANCER_COLOUR, // 4
    LABELLED_COLOUR, // 5
    APOPTOSIS_COLOUR, // 6
    INVISIBLE_COLOUR // visualizer treats '7' as invisible
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractCellPopulation)

//////////////////////////////////////////////////////////////////////////////
//         Iterator class implementation - most methods are inlined         //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPtr AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPtr AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator::operator!=(const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator& rOther)
{
    return mCellIter != rOther.mCellIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator& AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator::operator++()
{
    do
    {
        ++mCellIter;
    }
    while (!IsAtEnd() && !IsRealCell());

    return (*this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator::IsRealCell()
{
    return !( mrCellPopulation.IsCellAssociatedWithADeletedLocation(*mCellIter) || (*this)->IsDead() );
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrCellPopulation.rGetCells().end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator::Iterator(AbstractCellPopulation& rCellPopulation, std::list<CellPtr>::iterator cellIter)
    : mrCellPopulation(rCellPopulation),
      mCellIter(cellIter)
{
    // The cell population can now return empty if it only has ghost nodes
    if (mrCellPopulation.rGetCells().empty())
    {
        mCellIter = mrCellPopulation.rGetCells().end();
    }
    else
    {
        // Make sure we start at a real cell
        if (mCellIter == mrCellPopulation.rGetCells().begin() && !IsRealCell())
        {
            ++(*this);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Begin()
{
    return Iterator(*this, this->mCells.begin());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::End()
{
    return Iterator(*this, this->mCells.end());
}

#endif /*ABSTRACTCELLPOPULATION_HPP_*/
