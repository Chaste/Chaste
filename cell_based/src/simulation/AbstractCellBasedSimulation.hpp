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

#ifndef ABSTRACTCELLBASEDSIMULATION_HPP_
#define ABSTRACTCELLBASEDSIMULATION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <vector>

#include "AbstractCellKiller.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractForce.hpp"
#include "RandomNumberGenerator.hpp"

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCellPopulation;

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
        archive & mOutputDivisionLocations;
        archive & mOutputCellVelocities;
        archive & mCellKillers;
        archive & mSimulationModifiers;
        archive & mSamplingTimestepMultiple;
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

    /**
     * Whether to write the cell velocities to a file.
     * Initialised to false in constructor.
     */
    bool mOutputCellVelocities;

    /** Results file for cell velocities. */
    out_stream mpCellVelocitiesFile;

    /** List of cell killers. */
    std::vector<boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > > mCellKillers;

    /** List of SimulationModifier rules. */
    std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > > mSimulationModifiers;

    /**
     * The ratio of the number of actual timesteps to the number
     * of timesteps at which results are written to file.
     */
    unsigned mSamplingTimestepMultiple;

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
     * If mOutputDivisionLocations is set to true, then this method also writes the
     * location of each cell division at the present time to mpDivisionLocationFile.
     * This outputs a line of tab-separated values of the form:
     * [time] [div 0 x-pos] [div 0 y-pos] [div 0 z-pos] [div 0 age] [div 1 x-pos] [div 1 y-pos] [div 1 z-pos] [div 1 age] ...
     *
     * with [y-pos] and [z-pos] included for 2 and 3 dimensional simulations, respectively,
     * and [...age] denoting the age of the dividing cell.
     *
     * @return the number of births that occurred.
     */
    virtual unsigned DoCellBirth();

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
     * A method for subclasses to do something at before the start of the time loop.
     */
    virtual void SetupSolve()
    {
    }

    /**
     * A child class can overload this if they want the simulation to stop
     * based on certain conditions before the specified end time (for example,
     * run until a crypt becomes monoclonal).
     * @return true if stopping event has occurred
     */
    virtual bool StoppingEventHasOccurred();

    /**
     * Calls the deaths, births and (if mUpdateCellPopulation is true) CellPopulation::Update() methods.
     */
    virtual void UpdateCellPopulation();

    /**
     * Update the cell locations and topology (connectivity) of the cell population. This method
     * is called within the main time loop of Solve() .
     *
     * As this method is pure virtual, it must be overridden in subclasses.
     *
     * In the case of an OffLatticeSimulation, the method computes the force acting on each node
     * (corresponding to a cell in centre-based models and to a vertex in vertex-based models)
     * and integrates equations of motion to find the new position of each node.
     *
     * In the case of an OnLatticeSimulation, the method performs Monte Carlo updating of the cell
     * population, through a call to UpdateCellLocations() on the cell population object.
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
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputAdditionalSimulationSetup(out_stream& rParamsFile)=0;

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
     * @return the output directory of the simulation.
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
     * Method to remove all the cell killers.
     */
    void RemoveAllCellKillers();

    /**
     * Add a SimulationModifier to be used in this simulation.
     *
     * @param pSimulationModifier pointer to a SimulationModifier
     */
    void AddSimulationModifier(boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM> > pSimulationModifier);

    /**
     * @return a pointer to the vector of SimulationModifiers used in this simulation.
     */
    std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >* GetSimulationModifiers();

    /**
     * Main Solve() method, used to evolve the cell population. Note that prior to calling Solve()
     * we must have called SetEndTime(). We may also have optionally called SetDt(); if not, then
     * a default time step is used.
     *
     * The Solve() method proceeds as follows.
     *
     * Setting up:
     *
     * First, we set up SimulationTime, which (i) provides a globally consistent time, accessible
     * to all other classes in the cell_based code and (ii) handles any rounding issues when the
     * time step does not exactly divide the end time.
     *
     * Next, we create output files. We then call SetupSolve(), which is empty in the parent class
     * but may be overridden, e.g. to open additional output files. We then call SetupSolve() on
     * any member objects inheriting from AbstractCellBasedSimulationModifier. This class hierarchy
     * allows the user to introduce new ways of updating the cell population within the simulation.
     *
     * Next, we set up each cell by calling ReadyToDivide() on it, which updates the cell's age and
     * cell cycle model. Finally, we call WriteVisualizerSetupFile() and OutputSimulationSetup(),
     * as well as WriteResultsToFiles() on the cell population, to record the initial configration.
     * This completes the set up process.
     *
     * The main time loop:
     *
     * At each time step, we begin by calling UpdateCellPopulation(), which implements any cell
     * deaths and cell divisions through DoCellRemoval() and DoCellBirth() respectively. We then
     * update the correspondence between cells and the mesh by calling Update() on the cell
     * population.
     *
     * If mOutputCellVelocities is set to true and we are at a printing time, then we als write the
     * velocity of each cell at the present time to mpCellVelocitiesFile.
     * This outputs a line of tab-separated values of the form:
     * [time] [cell 0 x-pos] [cell 0 y-pos] [cell 0 z-pos] [cell 0 x-vel] [cell 0 y-vel] [cell 0 z-vel] ...
     *
     * with [y-pos] and [z-pos] included for 2 and 3 dimensional simulations, respectively,
     * and data for cells being ordered as given by the cell population Iterator.
     *
     * Next, we call UpdateCellLocationsAndTopology(), which is pure virtual in the parent class so
     * must be overridden. As the cell population has been updated, we then increment
     * SimulationTime by one time step. We then call UpdateAtEndOfTimeStep() on any
     * AbstractCellBasedSimulationModifiers present, e.g. to write additional output.
     * In an analogous manner to the calls to SetupSolve() prior to entering the main time loop.
     *
     * The last step within the main time loop is to output the present results to file.
     *
     * Finishing up:
     *
     * After exiting the main time loop, we call UpdateCellPopulation() in order to carry out a final
     * update of the cell population. We also call
     * UpdateAtEndOfSolve()} on any member objects inheriting from AbstractCellBasedSimulationModifier
     * in an analogous manner to the aforementioned calls to SetupSolve() UpdateAtEndOfTimeStep().
     * Finally, we close output files. This completes the Solve() method.
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
     * @return mOutputCellVelocities
     */
    bool GetOutputCellVelocities();

    /**
     * Set mOutputCellVelocities.
     *
     * @param outputCellVelocities the new value of mOutputCellVelocities
     */
    void SetOutputCellVelocities(bool outputCellVelocities);

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
