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

#ifndef SIMULATIONTIME_HPP_
#define SIMULATIONTIME_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>
#include "SerializableSingleton.hpp"
#include "TimeStepper.hpp"

/**
 * Simulation time object stores the simulation time.
 * It uses the singleton pattern to provide a globally consistent time.
 *
 * Note that the start time, end time and number of time steps must
 * be set before time can be incremented and returned.
 *
 * You should generally use the calls
 * IncrementTimeOneStep() and GetTime() when using this class.
 */
class SimulationTime : public SerializableSingleton<SimulationTime>
{
public:

    /**
     * @return a pointer to the simulation time object.
     * The first time this is called the simulation time object is created.
     */
    static SimulationTime* Instance();

    /**
     * Sets the end time and the number of time steps.
     * This must be called after SetStartTime() but before using any other methods.
     *
     * @param endTime time at which to end this run of the simulation
     * @param totalTimeStepsInSimulation  the number of time steps into which the above will be divided
     */
    void SetEndTimeAndNumberOfTimeSteps(double endTime, unsigned totalTimeStepsInSimulation);

    /**
     * Reset method for the end time and the number of time steps, to run the simulation
     * further after a first initial run.
     *
     * @param rEndTime the new end time for this simulation (the simulation will run from
     *      the current time to this new end time, NOT from 0 to this end time)
     * @param rNumberOfTimeStepsInThisRun the number of time steps into which the next run is split
     */
    void ResetEndTimeAndNumberOfTimeSteps(const double& rEndTime, const unsigned& rNumberOfTimeStepsInThisRun);

    /**
     * Get the simulation time step, set in earlier calls.
     *
     * Warning: Use of this method may result in round errors; generally use GetTime() instead.
     *
     * @return time step for this run of the simulation
     */
    double GetTimeStep() const;

    /**
     * Increment the simulation time by one time step.
     *
     * GetTime() will return an updated current time after this call.
     */
    void IncrementTimeOneStep();

    /**
     * Get the number of time steps that have elapsed.
     *
     * @return number of time steps which have been taken
     */
    unsigned GetTimeStepsElapsed() const;

    /**
     * Get the simulation time (in hours), should not have rounding errors.
     *
     * @return simulation time
     */
    double GetTime() const;

    /**
     * Destroy the current SimulationTime instance.  The next call to
     * Instance will create a new instance, on which
     * SetEndTimeAndNumberOfTimeSteps must be called again to reset time.
     *
     * This method *must* be called before program exit, to avoid a memory
     * leak.
     */
    static void Destroy();

    /**
     * Allows lower classes to check whether the simulation time class has been set up before using it
     *
     * @return whether the start time of the simulation has been set.
     */
    bool IsStartTimeSetUp() const;

    /**
     * Allows lower classes to check whether the simulation time class has been set up before using it
     *
     * @return whether the end time of the simulation and the number of timesteps has been set.
     */
    bool IsEndTimeAndNumberOfTimeStepsSetUp() const;

    /**
     * @return whether the simulation has finished.
     */
    bool IsFinished() const;

    /**
     * Set the start time of the simulation
     *
     * @param startTime the time at which the simulation begins (usually 0.0 hours)
     */
    void SetStartTime(double startTime);

protected:
    /**
     * Default simulation time constructor
     *
     * Sets up time, you must set the start time,
     * end time and number of time steps before using the object.
     */
    SimulationTime();

private:
    /**
     * A pointer to the singleton instance of this class.
     */
    static SimulationTime* mpInstance;

    /**
     * Delegate all time stepping to a TimeStepper class
     */
    static boost::shared_ptr<TimeStepper> mpTimeStepper;

    /**
     * Stores the time at which the simulation started
     */
    double mStartTime;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialization of a SimulationTime object must be done with care.
     * Do not serialize this singleton directly.  Instead, serialize
     * the object returned by GetSerializationWrapper.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mStartTime;
        archive & mpTimeStepper;
    }
};

#endif /*SIMULATIONTIME_HPP_*/
