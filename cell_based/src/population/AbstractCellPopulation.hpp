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
#include <boost/serialization/set.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <boost/foreach.hpp>

#include "AbstractMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "CellPropertyRegistry.hpp"
#include "Identifiable.hpp"
#include "AbstractCellPopulationCountWriter.hpp"
#include "AbstractCellPopulationWriter.hpp"
#include "AbstractCellWriter.hpp"

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCellBasedSimulation;

/**
 * An abstract facade class encapsulating a cell population.
 *
 * Contains a group of cells and associated methods.
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
        archive & mpCellPropertyRegistry;
        archive & mOutputResultsForChasteVisualizer;
        archive & mCellWriters;
        archive & mCellPopulationWriters;
        archive & mCellPopulationCountWriters;
    }

    /**
     * Open all files in mCellPopulationWriters and mCellWriters in append mode for writing.
     *
     * Files in mCellPopulationCountWriters are NOT opened in this call since they are not written in
     * a round-robin fashion.
     *
     * @param rOutputFileHandler handler for the directory in which to open this file.
     */
    void OpenRoundRobinWritersFilesForAppend(OutputFileHandler& rOutputFileHandler);

    /**
     * Close all files in mCellPopulationWriters and mCellWriters.
     *
     * Files in mCellPopulationCountWriters are NOT closed in this call since they are not written in
     * a round-robin fashion
     *
     */
    void CloseRoundRobinWritersFiles();

protected:

    /** Map location (node or VertexElement) indices back to cells. */
    std::map<unsigned, std::set<CellPtr> > mLocationCellMap;

    /** Map cells to location (node or VertexElement) indices. */
    std::map<Cell*, unsigned> mCellLocationMap;

    /** Reference to the mesh. */
    AbstractMesh<ELEMENT_DIM, SPACE_DIM>& mrMesh;

    /** List of cells. */
    std::list<CellPtr> mCells;

    /** Population centroid. */
    c_vector<double, SPACE_DIM> mCentroid;

    /** Meta results file for VTK. */
    out_stream mpVtkMetaFile;

    /** Cell property registry. */
    boost::shared_ptr<CellPropertyRegistry> mpCellPropertyRegistry;

    /** Whether to write results to file for visualization using the Chaste java visualizer (defaults to true). */
    bool mOutputResultsForChasteVisualizer;

    /** A list of cell writers. */
    std::vector<boost::shared_ptr<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> > > mCellWriters;

    /** A list of cell population writers. */
    std::vector<boost::shared_ptr<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> > > mCellPopulationWriters;

    /** A list of cell population count writers. */
    std::vector<boost::shared_ptr<AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM> > > mCellPopulationCountWriters;

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
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     */
    virtual void WriteVtkResultsToFile(const std::string& rDirectory)=0;

    /**
     * Constructor that just takes in a mesh.
     *
     * @param rMesh the mesh for the population.
     */
    AbstractCellPopulation(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     * Call #AcceptCellWriter across the whole population.
     *
     * By default the implementation here iterates over the cell population,
     * but this is overridden in some classes that need to go over nodes.
     */
    virtual void AcceptCellWritersAcrossPopulation();

public:

    /**
     * AbstractCellPopulation Constructor.
     *
     * @note Warning: the passed-in vector of cells will be emptied, even if the constructor
     * throws an exception!
     *
     * @param rMesh a reference to the mesh underlying the cell population
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
     * @param rDataName is the name associated with the data
     * @param dataValue is the value of the data, initially the same for each cell
     */
    void SetDataOnAllCells(const std::string& rDataName, double dataValue);

    /**
     * @return reference to the mesh, mrMesh.
     */
    AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rGetMesh();

    /**
     * @return a shared to a tetrahedral mesh, for use with a PDE modifier.
     * This method is called by AbstractGrowingDomainPdeModifier.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* GetTetrahedralMeshForPdeModifier()=0;

    /**
     * @param pdeNodeIndex index of a node in a tetrahedral mesh for use with a PDE modifier
     *
     * @return if a node, specified by its index in a tetrahedral mesh for use
     *         with a PDE modifier, is associated with a non-apoptotic cell.
     * This method can be called by PDE classes.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual bool IsPdeNodeAssociatedWithNonApoptoticCell(unsigned pdeNodeIndex);

    /**
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
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual double GetCellDataItemAtPdeNode(unsigned pdeNodeIndex,
                                            std::string& rVariableName,
                                            bool dirichletBoundaryConditionApplies=false,
                                            double dirichletBoundaryValue=0.0)=0;

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
     * Write any data necessary to a visualization setup file.
     * Used by AbstractCellBasedSimulation::WriteVisualizerSetupFile().
     *
     * @param pVizSetupFile a visualization setup file
     */
    virtual void WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile);

    /**
     * Add a new cell to the cell population.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pNewCell  the cell to add
     * @param pParentCell pointer to a parent cell (if required)
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way).
     */
    virtual CellPtr AddCell(CellPtr pNewCell, CellPtr pParentCell=CellPtr())=0;

    /**
     * @return a default value for the time step to use when simulating
     * the cell population.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * Note that the time step can be reset by calling SetDt() on the
     * simulation object used to simulate the cell population.
     */
    virtual double GetDefaultTimeStep()=0;

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
     * As this method is pure virtual, it must be overridden
     * in subclasses.
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
    std::vector<unsigned> GetCellCyclePhaseCount();

    /**
     * This counts the number of cells that the cell iterator covers. It does not include dead cells or cells that are
     * associated with a deleted location in the mesh.
     *
     * @return the number of real cells.
     */
    unsigned GetNumRealCells();

    /**
     * This returns the number of cells that are present in the internal mCells vector. It also includes dead cells and cells
     * that are associated with a deleted location in the mesh.
     *
     * @return the number of real cells.
     */
    unsigned GetNumAllCells();

    /**
     * Sets the Ancestor index of all the cells at this time to be the
     * same as their location index, can be used to trace clonal populations.
     */
    void SetCellAncestorsToLocationIndices();

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
    virtual CellPtr GetCellUsingLocationIndex(unsigned index);

    /**
     * Get the set of cells corresponding to a given location index.
     *
     * Note that the set may be empty.
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
     * @return whether there is a cell attached.
     */
    virtual bool IsCellAttachedToLocationIndex(unsigned index);

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
     * Given a cell, returns the set of location indices corresponding to neighbouring cells.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCell a cell
     * @return the set of neighbouring location indices.
     */
    virtual std::set<unsigned> GetNeighbouringLocationIndices(CellPtr pCell)=0;

    /**
     * @return the centroid of the cell population.
     */
    c_vector<double, SPACE_DIM> GetCentroidOfCellPopulation();

    /**
     * Update the ownership of cell in a parallel cell-based simulation.
     */
    virtual void UpdateCellProcessLocation();

    /**
     * Open output files (and, if required, write headers) for any writers in the members
     * mCellPopulationCountWriters, mCellPopulationWriters and mCellWriters.
     *
     * The method also writes the header for the .pvd output file if VTK is available.
     *
     * Before doing this, the method also creates appropriate writer objects if
     * mOutputResultsForChasteVisualizer is set to true.
     *
     * This method is public because it is called by the simulation class at the start of the Solve() call.
     *
     * @param rOutputFileHandler handler for the directory in which to open this file.
     */
    virtual void OpenWritersFiles(OutputFileHandler& rOutputFileHandler);

    /**
     * Close output files associated with any writers in the members
     * mCellPopulationCountWriters, mCellPopulationWriters and mCellWriters.
     *
     * The method also closes the .pvd output file if VTK is available.
     */
    void CloseWritersFiles();

    /**
     * Write results from the current cell population state to output files.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     */
    virtual void WriteResultsToFiles(const std::string& rDirectory);

    /**
     * Accept a cell population writer so it can write data from this object to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pPopulationWriter the population writer.
     */
    virtual void AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> > pPopulationWriter)=0;

    /**
     * Accept a cell population count writer so it can write data from this object to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pPopulationCountWriter the population count writer.
     */
    virtual void AcceptPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM> > pPopulationCountWriter)=0;

    /**
     * Accept a cell writer so it can write data from this object to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCellWriter the population writer.
     * @param pCell the cell whose data are being written.
     */
    virtual void AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> > pCellWriter, CellPtr pCell)=0;

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
     * Empty hook method to provide the ability to specify some additional property
     * of a cell-based simulation object.
     *
     * This method is called immediately prior to calling SetupSolve() within the
     * Solve() method in AbstractCellBasedSimulation.
     *
     * This method can be overridden, for example, to add a T2SwapCellKiller to the
     * simulation object in the case of a VertexBasedCellPopulation. This functionality
     * avoids the need for static or dynamic casts to specific cell population types
     * within simulation methods.
     *
     * @param pSimulation pointer to a cell-based simulation object
     */
    virtual void SimulationSetupHook(AbstractCellBasedSimulation<ELEMENT_DIM, SPACE_DIM>* pSimulation);

    /**
     * @return mOutputResultsForChasteVisualizer
     */
    bool GetOutputResultsForChasteVisualizer();

    /**
     * Add a cell population writer based on its type. Template parameters are inferred from the population.
     * The implementation of this function must be available in the header file.
     *
     * @return This method returns void
     */
    template<template <unsigned, unsigned> class T>
    void AddPopulationWriter()
    {
        mCellPopulationWriters.push_back(boost::shared_ptr< T<ELEMENT_DIM, SPACE_DIM> >(new T<ELEMENT_DIM, SPACE_DIM> ));
    }

    /**
     * Add a cell writer based on its type. Template parameters are inferred from the population.
     * The implementation of this function must be available in the header file.
     *
     * @return This method returns void
     */
    template<template <unsigned, unsigned> class T>
    void AddCellWriter()
    {
        mCellWriters.push_back(boost::shared_ptr< T<ELEMENT_DIM, SPACE_DIM> >(new T<ELEMENT_DIM, SPACE_DIM> ));
    }

    /**
     * Add a cell population count writer based on its type. Template parameters are inferred from the population.
     * The implementation of this function must be available in the header file.
     *
     * @return This method returns void
     */
    template<template <unsigned, unsigned> class T>
    void AddCellPopulationCountWriter()
    {
        mCellPopulationCountWriters.push_back(boost::shared_ptr< T<ELEMENT_DIM, SPACE_DIM> >(new T<ELEMENT_DIM, SPACE_DIM> ));
    }

    /**
     * Add a cell population writer through an input argument.
     * This alternative to the templated AddPopulationWriter()
     * method allows the user to, for example, add a writer
     * with a non-default value for its member mFileName.
     *
     * @param pPopulationWriter shared pointer to a cell population writer
     * @return This method returns void
     */
    void AddPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> > pPopulationWriter)
    {
        mCellPopulationWriters.push_back(pPopulationWriter);
    }

    /**
     * Add a cell writer through an input argument.
     * This alternative to the templated AddCellWriter()
     * method allows the user to, for example, add a writer
     * with a non-default value for its member mFileName.
     *
     * @param pCellWriter shared pointer to a cell writer
     * @return This method returns void
     */
    void AddCellWriter(boost::shared_ptr<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> > pCellWriter)
    {
        mCellWriters.push_back(pCellWriter);
    }

    /**
     * Add a cell population count writer through an input argument.
     * This alternative to the templated AddCellPopulationCountWriter()
     * method allows the user to, for example, add a writer
     * with a non-default value for its member mFileName.
     *
     * @param pCellPopulationCountWriter shared pointer to a cell population count writer
     * @return This method returns void
     */
    void AddCellPopulationCountWriter(boost::shared_ptr<AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM> > pCellPopulationCountWriter)
    {
        mCellPopulationCountWriters.push_back(pCellPopulationCountWriter);
    }

    /**
     * Get whether the population has a writer of the specified type.
     *
     * @return whether the population has this writer
     */
    template<template <unsigned, unsigned> class T>
    bool HasWriter() const
    {
        typedef AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> population_writer_t;
        BOOST_FOREACH(boost::shared_ptr<population_writer_t> p_pop_writer, mCellPopulationWriters)
        {
            if (dynamic_cast<T<ELEMENT_DIM, SPACE_DIM>* >(p_pop_writer.get()))
            {
                return true;
            }
        }
        typedef AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> cell_writer_t;
        BOOST_FOREACH(boost::shared_ptr<cell_writer_t> p_cell_writer, mCellWriters)
        {
            if (dynamic_cast<T<ELEMENT_DIM, SPACE_DIM>* >(p_cell_writer.get()))
            {
                return true;
            }
        }
        typedef AbstractCellPopulationCountWriter<ELEMENT_DIM, SPACE_DIM> count_writer_t;
        BOOST_FOREACH(boost::shared_ptr<count_writer_t> p_count_writer, mCellPopulationCountWriters)
        {
            if (dynamic_cast<T<ELEMENT_DIM, SPACE_DIM>* >(p_count_writer.get()))
            {
                return true;
            }
        }
        return false;
    }

    /**
     * Set mOutputResultsForChasteVisualizer.
     *
     * @param outputResultsForChasteVisualizer the new value of mOutputResultsForChasteVisualizer
     */
    void SetOutputResultsForChasteVisualizer(bool outputResultsForChasteVisualizer);

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
     * @return a pair of indices ordered by node index.
     * Used by the rest length routines.
     *
     * @param index1 a node index
     * @param index2 a node index
     */
    std::pair<unsigned,unsigned> CreateOrderedPair(unsigned index1, unsigned index2);

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
        inline bool operator!=(const typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator& rOther);

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
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator::operator!=(const typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator& rOther)
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
