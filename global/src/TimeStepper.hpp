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

#ifndef TIMESTEPPER_HPP_
#define TIMESTEPPER_HPP_

#include <vector>
#include <climits>
#include <boost/serialization/vector.hpp>

#include "ChasteSerialization.hpp"
#include "Exception.hpp"

/**
 * A helper class that provides a robust way to deal with time loops.
 *
 * An incremented integer counter is used to calculate the current time
 * and ensure that the correct end time.
 */
class TimeStepper
{
    friend class TestTimeStepper;
private:
    /*
     * Private default constructor for archiving
     */
     TimeStepper(){};

public:

    /**
     * Create a new time stepper over some simulation interval.
     * Time units are not specified, but all parameters should have consistent units.
     *
     * @param startTime  the start of the interval
     * @param endTime  the end of the interval
     * @param dt  the time step to use.
     * @param enforceConstantTimeStep If this is true the stepper estimates whether non-constant
     *  timesteps will be used and quits if so.
     * @param additionalTimes If the timestepper needs to stop at certain additional times, they can be passed in in this std::vector.
     *  Defaults to empty. These times must be in ascending order.
     *  DEPRECATION: Note that these additional times are checked in order to ensure that the stepper will
     *  stop at these times anyway.  (For example we want to check that new events such as electrodes switching on will be detected
     *  in a regular time loop.) New additional times are not added but instead throw an exception.
     */
    TimeStepper(double startTime, double endTime, double dt, bool enforceConstantTimeStep=false, std::vector<double> additionalTimes=std::vector<double> ());

    /**
     * Step forward one step in time and update member variables.
     *
     * Note that this should be used in conjuction with (not) IsTimeAtEnd because
     * increment past the end of a simulation will throw an exception
     */
    void AdvanceOneTimeStep();

    /**
     * @return Get the time.
     */
    double GetTime() const;

    /**
     * @return Get the value time will take at the next step.
     */
    double GetNextTime() const;

    /**
     * @return Get the size of the next time step which will be taken.
     *
     * Note that this is often reported as the ideal timestep.
     *
     * The actual time step is mNextTime - mTime which could be evaluated as
     *  GetNextTime() - GetTime() = mNextTime - mTime
     *  ~= mStart + (mTotalTimeStepsTaken  + 1)*mDt - mTime
     *  ~= mStart + (mTotalTimeStepsTaken  + 1)*mDt - (mStart + mTotalTimeStepsTaken*mDt)
     *  ~= mDt
     *  This wouldn't evaluate to mDt in general but to within 2*mEpsilon <= 2*DBL_EPSILON*mEnd*mDt
     *  When mEnd and/or mStart are a long way from zero then the value reported
     *   mNextTime - mTime
     *  will differ from mDt, even when the actual step is mDt
     */
    double GetNextTimeStep();

    /**
     * @return Get the size of the ideal time step (as set in the constructor.
     * Note when a fixed time step is used GetNextTimeStep() ==  GetIdealTimeStep()
     * until the end of the simulation when GetNextTimeStep() == 0.0
     */
    double GetIdealTimeStep();

    /**
     * @return True when GetTime == endTime.
     */
    bool IsTimeAtEnd() const;

    /**
     * @return Estimate number of time steps remaining, which may be an overestimate.
     * Used to reserve memory for writing intermediate solutions.
     */
    unsigned EstimateTimeSteps() const;

    /**
     * @return The number of times AdvanceOneTimeStep() has been called SINCE
     * the last time ResetTimeStep() was called.
     */
    unsigned GetTotalTimeStepsTaken() const;

    /**
     * Set the time step to use for adaptive time integration. Note that
     * this also resets mStart to be the current time and zeroes
     * mTotalTimeStepsTaken.
     *
     * @param dt  is the new time-step
     */
    void ResetTimeStep(double dt);

private:

    /** The start time. */
    double mStart;

    /** The end time. */
    double mEnd;

    /** The size of time step. */
    double mDt;

    /** The total number of time steps taken. */
    unsigned mTotalTimeStepsTaken;

    /** DEPRECATED: The number of times one of the 'additional times' has been passed. */
    unsigned mAdditionalTimesReachedDeprecated;

    /** The current time. */
    double mTime;

    /** What the time will be after the next time step. */
    double mNextTime;

    /** An architecture-dependent scaling factor.  This is so that we can compare
     * relative to the end time when mEnd is large.
     * mEpsilon = DBL_EPSILON when mEnd is small
     *          = mEnd*DBL_EPSILON when mEnd is large
     */
    double mEpsilon;

    /** @return Compute what the time will be at the next time step. */
    double CalculateNextTime();

    /** DEPRECATED Vector to store the additional times the stepper must stop at. */
    std::vector<double> mAdditionalTimesDeprecated;

        /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mStart;
        archive & mEnd;
        archive & mDt;
        archive & mTotalTimeStepsTaken;
        archive & mTime;
        archive & mNextTime;
        archive & mEpsilon;
        archive & mAdditionalTimesReachedDeprecated;
        archive & mAdditionalTimesDeprecated;
    }
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TimeStepper)

#endif /*TIMESTEPPER_HPP_*/
