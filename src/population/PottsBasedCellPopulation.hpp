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

#ifndef POTTSBASEDCELLPOPULATION_HPP_
#define POTTSBASEDCELLPOPULATION_HPP_

#include "AbstractOnLatticeCellPopulation.hpp"
#include "PottsMesh.hpp"
#include "VertexMesh.hpp"
#include "AbstractPottsUpdateRule.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"

template<unsigned DIM>
class AbstractPottsUpdateRule; // Circular definition

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

    /** Potts-based mesh associated with the cell population. */
    PottsMesh<DIM>& mrMesh;

    /**
     * Pointer to a VertexMesh object that stores the Element tessellation that is used to
     * visualise mrMesh. The tessellation is created by calling CreateElementTessellation()
     * and can be accessed by calling GetElementTessellation().
     */
    VertexMesh<DIM,DIM>* mpElementTessellation;

    /** Results file for elements. */
    out_stream mpVizElementsFile;

    /** The update rules used to determine the new location of the cells. */
    std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > > mUpdateRuleCollection;

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
#define COVERAGE_IGNORE
        archive & boost::serialization::base_object<AbstractOnLatticeCellPopulation<DIM> >(*this);

        /*
         * In its current form the code does not allow the direct serialization
         * of the VertexMesh class, so instead we delete mpVoronoiTessellation.
         */
        delete mpElementTessellation;
        mpElementTessellation = NULL;

        archive & mUpdateRuleCollection;
        archive & mTemperature;
        archive & mNumSweepsPerTimestep;

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
     * @param rCellDivisionVector  if this vector has any non-zero component, then it is used as the axis
     *     along which the parent cell divides
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
     * Overridden GenerateCellResultsAndWriteToFiles() method.
     */
    virtual void GenerateCellResultsAndWriteToFiles();

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
     * Add an update rule to be used in this simulation (use this to set up the Hamiltonian).
     *
     * @param pUpdateRule pointer to an update rule
     */
    void AddUpdateRule(boost::shared_ptr<AbstractPottsUpdateRule<DIM> > pUpdateRule);

    /**
     * Get the collection of update rules to be used in this simulation.
     *
     * @return the update rule collection
     */
    const std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > >& rGetUpdateRuleCollection() const;

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
     * Get a reference to mpElementTessellation.
     */
    VertexMesh<DIM,DIM>* GetElementTessellation();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsBasedCellPopulation)

// No archiving yet so untested
#define COVERAGE_IGNORE
namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PottsBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const PottsBasedCellPopulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
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
#undef COVERAGE_IGNORE

#endif /*POTTSBASEDCELLPOPULATION_HPP_*/
