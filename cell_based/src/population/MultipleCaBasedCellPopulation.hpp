/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef MULTIPLECABASEDCELLPOPULATION_HPP_
#define MULTIPLECABASEDCELLPOPULATION_HPP_

#include "AbstractOnLatticeCellPopulation.hpp"
#include "PottsMesh.hpp"
#include "VertexMesh.hpp"
#include "AbstractMultipleCaUpdateRule.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"

///\todo Is this needed in this case?
template<unsigned DIM>
class AbstractMultipleCaUpdateRule; // Circular definition

/**
 * \todo This description is pasted from Potts and needs to be changed (#2066)
 *
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
class MultipleCaBasedCellPopulation : public AbstractOnLatticeCellPopulation<DIM>
{
    friend class TestMultipleCaBasedCellPopulation;

private:

    /** The carrying capacity (number of cells allowed per site). */
    unsigned mLatticeCarryingCapacity;

    /** Results file for cell locations. */
    out_stream mpVizLocationsFile;

    /** The update rules used to determine the new location of the cells. */
    std::vector<boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > > mUpdateRuleCollection;

    /** Records for each node the node the number of spaces available. */
    std::vector<unsigned> mAvailableSpaces;

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
#define COVERAGE_IGNORE
        archive & boost::serialization::base_object<AbstractOnLatticeCellPopulation<DIM> >(*this);
        archive & mLatticeCarryingCapacity;
        archive & mUpdateRuleCollection;
        archive & mAvailableSpaces;
#undef COVERAGE_IGNORE
    }

    /**
     * Check the consistency of internal data structures.
     * Each PottsElement must have a CellPtr associated with it.
     */
    void Validate();

    /**
     * Overridden WriteVtkResultsToFile() method.
     */
    void WriteVtkResultsToFile();

public:

    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely one CellPtr for each PottsElement in
     * the mesh.
     *
     * @param rMesh reference to a PottsMesh
     * @param rCells reference to a vector of CellPtrs
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param latticeCarryingCapacity an optional parameter to allow more than one cell per site
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     *                   (defaults to false)
     * @param validate whether to validate the cell population when it is created (defaults to false as not used in CA simulations)
     */
    MultipleCaBasedCellPopulation(PottsMesh<DIM>& rMesh,
                                  std::vector<CellPtr>& rCells,
                                  const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                                  unsigned latticeCarryingCapacity=1u,
                                  bool deleteMesh=false,
                                  bool validate=false);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a vertex mesh.
     */
    MultipleCaBasedCellPopulation(PottsMesh<DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~MultipleCaBasedCellPopulation();

    /**
     * @return mAvailableSpaces.
     */
    std::vector<unsigned>& rGetAvailableSpaces();

    /**
     * Find if a given node Has space available.
     *
     * @param index the global index of a specified node
     *
     * @return whether the node is an empty site
     */
    bool IsSiteAvailable(unsigned index);

    /**
     * @return reference to #mrMesh.
     */
    PottsMesh<DIM>& rGetMesh();

    /**
     * @return const reference to #mrMesh (used in archiving).
     */
    const PottsMesh<DIM>& rGetMesh() const;

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
     * Overridden GetLocationOfCellCentre() method.
     * Find where a given cell is in space.
     *
     * @param pCell the cell
     *
     * @return the location of the centre of mass of the element corresponding to this cell.
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
     * @param rCellDivisionVector  this parameter is not yet used in this class (see #1737)
     * @param pParentCell pointer to a parent cell (if required)
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    CellPtr AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell=CellPtr());

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
     * Overridden WriteCellVolumeResultsToFile() method.
     */
    void WriteCellVolumeResultsToFile();

    /**
     * Overridden GetVolumeOfCell() method.
     *
     * @param pCell boost shared pointer to a cell
     */
    double GetVolumeOfCell(CellPtr pCell);

    /**
     * Overridden GenerateCellResultsAndWriteToFiles() method.
     */
    virtual void GenerateCellResultsAndWriteToFiles();

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
     * Add an update rule to be used in this simulation (use this to set up how cells move).
     *
     * @param pUpdateRule pointer to an update rule
     */
    void AddUpdateRule(boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > pUpdateRule);

    /**
     * Method to remove all the update rules
     */
    void RemoveAllUpdateRules();

    /**
     * Get the collection of update rules to be used in this simulation.
     *
     * @return the update rule collection
     */
    const std::vector<boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > >& rGetUpdateRuleCollection() const;

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
     * Overridden GetNeighbouringNodeIndices() method.
     *
     * This method currently returns an exception as the two types of neighbourhood
     * (Moore and Von Neumann) are defined in the PottsMesh.
     *
     * @param index the node index
     * @return the set of neighbouring node indices.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MultipleCaBasedCellPopulation)

// No archiving yet so untested
#define COVERAGE_IGNORE
namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MultipleCaBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MultipleCaBasedCellPopulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const PottsMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a MultipleCaBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MultipleCaBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    PottsMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)MultipleCaBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...
#undef COVERAGE_IGNORE

#endif /*MULTIPLECABASEDCELLPOPULATION_HPP_*/
