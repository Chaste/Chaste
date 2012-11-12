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

#ifndef CELLBASEDPDEHANDLER_HPP_
#define CELLBASEDPDEHANDLER_HPP_

#include <map>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>

#include "AbstractCellPopulation.hpp"
#include "PdeAndBoundaryConditions.hpp"
#include "TetrahedralMesh.hpp"
#include "ChasteCuboid.hpp"
#include "Identifiable.hpp"

/**
 * A helper class, containing code for handling the numerical solution of one or more PDEs
 * (using the finite element method) associated with a cell-based simulation object.
 *
 * By letting AbstractCellBasedSimulation have a pointer to an object of this type as a
 * member variable, we separate out all PDE-related functionality into this class, and thus
 * obviate the need for specialized cell-based simulation subclasses.
 */
template<unsigned DIM>
class CellBasedPdeHandler : public Identifiable
{
    // Allow tests to access private members, in order to test computation of private functions
    friend class TestCellBasedPdeHandler;
    friend class TestOffLatticeSimulationWithPdes;
    friend class TestOnLatticeSimulationWithPdes;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mPdeAndBcCollection;
        archive & mWriteAverageRadialPdeSolution;
        archive & mWriteDailyAverageRadialPdeSolution;
        archive & mSetBcsOnCoarseBoundary;
        archive & mNumRadialIntervals;
        archive & mAverageRadialSolutionVariableName;
    }

protected:

    /** Pointer to a cell population. */
    AbstractCellPopulation<DIM>* mpCellPopulation;

    /** Vector of pointers to linear elliptic PDE objects with additional boundary condition information. */
    std::vector<PdeAndBoundaryConditions<DIM>*> mPdeAndBcCollection;

    /** A cache of where the results are going (used for VTK writer). */
    std::string mDirPath;

    /** File that the values of the PDE solution are written out to. */
    out_stream mpVizPdeSolutionResultsFile;

    /** File that the average radial PDE solution is written out to. */
    out_stream mpAverageRadialPdeSolutionResultsFile;

    /** Whether to write to file the average radial PDE solution. */
    bool mWriteAverageRadialPdeSolution;

    /** Whether to write the average radial PDE solution daily. */
    bool mWriteDailyAverageRadialPdeSolution;

    /** The name of the quantity that gets averaged. */
    std::string mAverageRadialSolutionVariableName;

    /** Whether to set the boundary condition on the edge of the coarse mesh rather than the cell population. */
    bool mSetBcsOnCoarseBoundary;

    /** Number of radial 'bins' used to calculate the average radial PDE solution. */
    unsigned mNumRadialIntervals;

    /** Coarse mesh on which to solve the PDE. */
    TetrahedralMesh<DIM,DIM>* mpCoarsePdeMesh;

    /** Map between cells and the elements of the coarse PDE mesh containing them. */
    std::map<CellPtr, unsigned> mCellPdeElementMap;

    /**
     * Whether to delete member pointers in the destructor.
     * Used in archiving.
     */
    bool mDeleteMemberPointersInDestructor;

    /**
     * Initialise mCellPdeElementMap.
     *
     * This method is only called Solve(), but is written as a separate method
     * for testing purposes.
     */
    void InitialiseCellPdeElementMap();

    /**
     * Write the PDE solution to file at a specified time.
     *
     * @param time The time at which to record the PDE solution
     */
    virtual void WritePdeSolution(double time);

    /**
     * Write the average radial PDE solution to file at a specified time.
     *
     * @param time The time at which to record the average radial PDE solution
     */
    void WriteAverageRadialPdeSolution(double time);

    /**
     * @return Whether the population being used requires the use of separate `coarse` mesh to
     * solve PDEs
     */
    bool PdeSolveNeedsCoarseMesh();

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to a cell population
     * @param deleteMemberPointersInDestructor whether to delete member pointers in the destructor (defaults to false)
     */
    CellBasedPdeHandler(AbstractCellPopulation<DIM>* pCellPopulation, bool deleteMemberPointersInDestructor=false);

    /**
     * Destructor.
     */
    virtual ~CellBasedPdeHandler();

    /**
     * Get a pointer to the cell population.
     *
     * @return a const pointer to mpCellPopulation
     */
    const AbstractCellPopulation<DIM>* GetCellPopulation() const;

    /**
     * @return mpCoarsePdeMesh
     */
    TetrahedralMesh<DIM,DIM>* GetCoarsePdeMesh();

    /**
     * Open results files and write initial conditions to file.
     * Called by AbstractCellBasedSimulation::Solve().
     *
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    void OpenResultsFiles(std::string outputDirectory);

    /**
     * Close results files.
     * Called by AbstractCellBasedSimulation::Solve().
     */
    void CloseResultsFiles();

    /**
     * @return mWriteAverageRadialPdeSolution
     */
    bool GetWriteAverageRadialPdeSolution();

    /**
     * @return mWriteDailyAverageRadialPdeSolution
     */
    bool GetWriteDailyAverageRadialPdeSolution();

    /**
     * Update the mCellPdeElementMap
     *
     * This method should be called before sending the element map to a PDE class
     * to ensure map is up to date.
     */
    void UpdateCellPdeElementMap();

    /**
     * @return mSetBcsOnCoarseBoundary
     */
    bool GetImposeBcsOnCoarseBoundary();

    /**
     * @return mNumRadialIntervals
     */
    unsigned GetNumRadialIntervals();

    /**
     * Solve the PDE and write the solution to file.
     *
     * @param samplingTimestepMultiple the ratio of the number of actual timesteps to the number of timesteps
     *     at which results are written to file.
     */
    virtual void SolvePdeAndWriteResultsToFile(unsigned samplingTimestepMultiple);

    /**
	 * Find the solution of one of the PDEs at a point in space
	 *
	 * @param Point the position in space
	 * @param Variable the dependent variable of the PDE whose solution you want to find
	 *
	 * @return the solution of the required PDE at the given point.
	 */
    double GetPdeSolutionAtPoint(c_vector<double,DIM> Point, std::string Variable);

    /**
     * Find the index of the coarse mesh element containing a given cell.
     *
     * @param pCell the cell
     *
     * @return the element index.
     */
    unsigned FindCoarseElementContainingCell(CellPtr pCell);

    /**
     * Get the solution to the PDE at this time step.
     *
     * @param rName The name of the dependent variable for the PDE in the vector mPdeAndBcCollection.
     * This defaults to an empty string in the case there is only one PDE.
     */
    virtual Vec GetPdeSolution(const std::string& rName = "");

    /**
     * Write the final (and optionally also the daily) average
     * radial PDE solution to file.
     *
     * @param rName The name of the quantity that we are averaging.
     * @param numRadialIntervals The number of radial intervals in which the average
     *                           PDE solution is calculated (defaults to 10)
     * @param writeDailyResults Whether to record the average radial PDE solution
     *                          at the end of each day of the simulation (defaults to false)
     */
    void SetWriteAverageRadialPdeSolution(const std::string& rName,
                                            unsigned numRadialIntervals=10,
                                            bool writeDailyResults=false);

    /**
     * Impose the PDE boundary conditions on the edge of the cell population when using
     * the coarse mesh. The default option is to impose the condition on the boundary of the
     * coarse mesh.
     *
     * @param setBcsOnCoarseBoundary whether to impose the BCs on the edge of the cell population
     */
    void SetImposeBcsOnCoarseBoundary(bool setBcsOnCoarseBoundary);

    /**
     * Solve the PDE problem on a coarse mesh.
     *
     * @param stepSize horizontal and vertical distance between mesh points
     * @param meshCuboid the cuboid defining the size and location of the mesh.
     * @param centreOnCellPopulation whether to centre the coarse mesh on the centre of the cell population.
     */
    virtual void UseCoarsePdeMesh(double stepSize, ChasteCuboid<DIM> meshCuboid, bool centreOnCellPopulation = false);

    /**
     * Pass a PDE and associated boundary conditions to the simulation.
     *
     * @param pPdeAndBc a pointer to a PdeAndBoundaryConditions object
     */
    void AddPdeAndBc(PdeAndBoundaryConditions<DIM>* pPdeAndBc);

    /**
     * Output parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellBasedPdeHandler)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CellBasedPdeHandler.
 *
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CellBasedPdeHandler<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = t->GetCellPopulation();
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a CellBasedPdeHandler.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CellBasedPdeHandler<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CellBasedPdeHandler<DIM>(p_cell_population, true);
}
}
}

#endif /*CELLBASEDPDEHANDLER_HPP_*/
