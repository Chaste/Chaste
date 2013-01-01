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

#ifndef ABSTRACTCELLBASEDSIMULATION_HPP_
#define ABSTRACTCELLBASEDSIMULATION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <vector>

#include "CellBasedPdeHandler.hpp"
#include "AbstractCellKiller.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * An abstract cell-based simulation class. This class contains common functionality
 * from off lattice and on lattice simulations.
 *
 * The AbstractCellBasedSimulation is constructed with a CellPopulation, which
 * updates the correspondence between each Cell and its spatial representation
 * and handles cell division (governed by the CellCycleModel associated
 * with each cell). Once constructed,  one or more CellKillers may be passed
 * to the AbstractCellBasedSimulation object to specify conditions in which Cells
 * may die,
 *
 * Subclasses use one or more Force laws or update rules (Which are passed
 * to the child class object) to define the mechanical properties of the CellPopulation.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractCellBasedSimulation : public Identifiable
{
    // Allow tests to access private members, in order to test computation of private functions e.g. DoCellBirth()
    friend class TestCryptSimulation2d;
    friend class TestOffLatticeSimulation3d;
    friend class TestOffLatticeSimulation;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Save or restore the simulation.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        SerializableSingleton<SimulationTime>* p_time_wrapper = SimulationTime::Instance()->GetSerializationWrapper();
        archive & p_time_wrapper;

        SerializableSingleton<RandomNumberGenerator>* p_rng_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_rng_wrapper;

        archive & mDt;
        archive & mEndTime;
        archive & mNoBirth;
        archive & mUpdateCellPopulation;
        archive & mOutputDirectory;
        archive & mNumBirths;
        archive & mNumDeaths;
        archive & mCellKillers;
        archive & mSamplingTimestepMultiple;
        archive & mpCellBasedPdeHandler;
    }

protected:

    /** Time step. */
    double mDt;

    /** Time to run the Solve() method up to. */
    double mEndTime;

    /** Facade encapsulating cells in the cell population being simulated. */
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& mrCellPopulation;

    /** Whether to delete the cell population in the destructor. */
    bool mDeleteCellPopulationInDestructor;

    /** Whether to initialise the cells. */
    bool mInitialiseCells;

    /** Whether to run the simulation with no birth (defaults to false). */
    bool mNoBirth;

    /** Whether to update the topology of the cell population at each time step (defaults to true).*/
    bool mUpdateCellPopulation;

    /** Output directory (a subfolder of tmp/[USERNAME]/testoutput). */
    std::string mOutputDirectory;

    /** Simulation Output directory either the same as mOutputDirectory or includes mOutputDirectory/results_from_time_[TIME]. */
    std::string mSimulationOutputDirectory;

    /** Visualizer setup file. */
    out_stream mpVizSetupFile;

    /** Counts the number of births during the simulation. */
    unsigned mNumBirths;

    /** Counts the number of deaths during the simulation. */
    unsigned mNumDeaths;

    /** Whether to output the locations of division events (defaults to false).*/
    bool mOutputDivisionLocations;

    /** Output file for location of division events. */
    out_stream mpDivisionLocationFile;

    /** List of cell killers. */
    std::vector<boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > > mCellKillers;

    /**
     * The ratio of the number of actual timesteps to the number
     * of timesteps at which results are written to file.
     */
    unsigned mSamplingTimestepMultiple;

    /**
     * Pointer to a CellBasedPdeHandler object.
     */
    CellBasedPdeHandler<SPACE_DIM>* mpCellBasedPdeHandler;

    /**
     * Writes out special information about the mesh to the visualizer.
     */
    virtual void WriteVisualizerSetupFile()
    {
    }

    /**
     * During a simulation time step, process any cell divisions that need to occur.
     * If the simulation includes cell birth, causes (almost) all cells that are ready to divide
     * to produce daughter cells.
     *
     * @return the number of births that occurred.
     */
    virtual unsigned DoCellBirth();

    /**
     * Method for determining how cell division occurs. This method returns a vector
     * which is then passed into the CellPopulation method AddCell().
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pParentCell the parent cell
     *
     * @return a vector containing information on cell division.
     */
    virtual c_vector<double, SPACE_DIM> CalculateCellDivisionVector(CellPtr pParentCell)=0;

    /**
     * During a simulation time step, process any cell sloughing or death
     *
     * This uses the cell killers to remove cells and associated nodes from the
     * facade class.
     *
     * @return the number of deaths that occurred.
     */
    unsigned DoCellRemoval();

    /**
     * A method for subclasses to do something at the end of each timestep.
     */
    virtual void UpdateAtEndOfTimeStep()
    {
    }

    /**
     * A method for subclasses to do something at before the start of the time loop.
     */
    virtual void SetupSolve()
    {
    }

    /**
     * This method may be overridden in subclasses to do something
     * at the end of each time loop.
     */
    virtual void UpdateAtEndOfSolve()
    {
    }

    /**
     * A child class can overload this if they want the simulation to stop
     * based on certain conditions before the specified end time (for example,
     * run until a crypt becomes monoclonal).
     */
    virtual bool StoppingEventHasOccurred();

    /**
     * Calls the deaths, births and (if mUpdateCellPopulation is true) CellPopulation::Update() methods.
     */
    virtual void UpdateCellPopulation();

    /**
     * Updates the locations and topology of cells, must be defined in subclasses.
     *
     * For example this calculates forces and updates node positions in off lattice simulations.
     */
    virtual void UpdateCellLocationsAndTopology()=0;

    /**
     * Helper method to output all the simulations parameters and information to file.
     */
    void OutputSimulationSetup();

    /**
     * Helper method to output additional simulations parameters and information defined in
     * subclasses to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputAdditionalSimulationSetup(out_stream& rParamsFile)
    {
    }

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells Whether to initialise cells (defaults to true; set to false when loading from an archive)
     */
    AbstractCellBasedSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                bool deleteCellPopulationInDestructor=false,
                                bool initialiseCells=true);

    /**
     * Destructor.
     *
     * This frees the cell population if it was created by de-serialization.
     */
    virtual ~AbstractCellBasedSimulation();

    /**
     * Set mpCellBasedPdeHandler
     *
     * @param pCellBasedPdeHandler pointer to a CellBasedPdeHandler object
     */
    void SetCellBasedPdeHandler(CellBasedPdeHandler<SPACE_DIM>* pCellBasedPdeHandler);

    /**
     * @return mpCellBasedPdeHandler
     */
    CellBasedPdeHandler<SPACE_DIM>* GetCellBasedPdeHandler();

    /**
     * Get a node's location (ONLY FOR TESTING).
     *
     * @param rNodeIndex the node index
     * @return the co-ordinates of this node.
     */
    std::vector<double> GetNodeLocation(const unsigned& rNodeIndex);

    /**
     * @return the timestep of the simulation
     */
    double GetDt();

    /**
     * @return the number of births that have occurred in the entire simulation (since t=0)
     */
    unsigned GetNumBirths();

    /**
     * @return the number of deaths that have occurred in the entire simulation (since t=0).
     */
    unsigned GetNumDeaths();

    /**
     * Get the output directory of the simulation.
     */
    std::string GetOutputDirectory();

    /**
     * Set the timestep of the simulation.
     *
     * @param dt the timestep to use
     */
    void SetDt(double dt);

    /**
     * Set the end time and resets the timestep to be endtime/100.
     *
     * @param endTime the end time to use
     */
    void SetEndTime(double endTime);

    /**
     * Set the output directory of the simulation.
     *
     * @param outputDirectory the output directory to use
     */
    void SetOutputDirectory(std::string outputDirectory);

    /**
     * Set the ratio of the number of actual timesteps to the number of timesteps
     * at which results are written to file. Default value is set to 1 by the constructor.
     *
     * @param samplingTimestepMultiple the ratio to use
     */
    void SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple);

    /**
     * Set the simulation to run with no birth.
     *
     * @param noBirth whether to run with no birth
     */
    void SetNoBirth(bool noBirth);

    /**
     * Set whether to update the topology of the cell population at each time step.
     *
     * @param updateCellPopulation  whether to update the cell population each time step
     */
    void SetUpdateCellPopulationRule(bool updateCellPopulation);

    /**
     * Return whether to update the topology of the cell population at each time step.
     *
     * @return whether to update the cell population each time step
     */
    bool GetUpdateCellPopulationRule();

    /**
     * Add a cell killer to be used in this simulation.
     *
     * @param pCellKiller pointer to a cell killer
     */
    void AddCellKiller(boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > pCellKiller);

    /**
     * Method to remove all the cell killers
     */
    void RemoveAllCellKillers();

    /**
     * Main solve method.
     *
     * This method sets up the simulation time, creates output files, and initialises the
     * cell population. It then iterates through a time loop. At each time step, first any cell death
     * or birth is implemented, then the cell population topology is updated, then the forces are
     * recalculated and the cell population evolved according to whatever force laws are present in
     * the simulation, and finally the results for that time step are output to file. At the
     * end of the time loop, the method closes any output files.
     */
    void Solve();

    /**
     * @return reference to the cell population.
     */
    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rGetCellPopulation();

    /**
     * @return const reference to the cell population (used in archiving).
     */
    const AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rGetCellPopulation() const;

    /**
     * @return mOutputDivisionLocations
     */
    bool GetOutputDivisionLocations();

    /**
     * Set mOutputDivisionLocations.
     *
     * @param outputDivisionLocations the new value of mOutputDivisionLocations
     */
    void SetOutputDivisionLocations(bool outputDivisionLocations);

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationParameters(out_stream& rParamsFile)=0;
};

#endif /*ABSTRACTCELLBASEDSIMULATION_HPP_*/
