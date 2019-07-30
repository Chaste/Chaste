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

#ifndef _CELLPROPERTIES_HPP_
#define _CELLPROPERTIES_HPP_

#include <vector>

/**
 * Class to calculate various physiological properties from the results of a
 * cardiac simulation.
 *
 * It will calculate for a single cell:
 *   All the Action potential durations (at any percentage)
 *   Max. upstroke velocity for each AP
 *   Action potential amplitudes for each AP
 *   Peak & Resting membrane potentials for all AP
 *   Cycle lengths (time between onset of APs)
 *
 */
class CellProperties
{
private:
    /**  The simulation results to process, voltage */
    const std::vector<double>& mrVoltage;
    /**  The simulation results to process, time */
    const std::vector<double>& mrTime;

    /**
     * Whether we have any unfinished action potentials
     * i.e. the trace ends ABOVE_THRESHOLD.
     */
    bool mUnfinishedActionPotentials;

    /**
     * Threshold for determining what counts as an action potential.
     * This is a value part way between the min & max potential, to avoid
     * problems due to 'notches' in an action potential.
     */
    double mThreshold;

    /** Cached vector containing AP onset properties */
    std::vector<double> mOnsets;

    /** Cached vector containing AP resting properties */
    std::vector<double> mRestingValues;

    /** Cached vector containing AP cycle lengths properties */
    std::vector<double> mCycleLengths;

    /** Cached vector containing AP peak properties */
    std::vector<double> mPeakValues;

    /** Cached vector containing the times of AP peaks */
    std::vector<double> mTimesAtPeakValues;

    /** Cached vector containing AP upstroke properties */
    std::vector<double> mMaxUpstrokeVelocities;

    /** Cached vector containing the times of AP upstrokes */
    std::vector<double> mTimesAtMaxUpstrokeVelocity;

    /** Cached vector containing the number of recorded depolarisations while above threshold */
    std::vector<unsigned> mCounterOfPlateauDepolarisations;

    /**
     * Calculate all the cacheable values.
     */
    void CalculateProperties();

    /**
     * Actually calculate APD.
     *
     * target voltage = resting +'percentage'*(amplitude-resting).
     *
     * APD is taken to be the time from when the target voltage is crossed upwards to the time
     * voltage crosses the target voltage downwards.
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     *
     * @return a vector containing the APDs.
     */
    std::vector<double> CalculateActionPotentialDurations(const double percentage);

    /**
     * Throw an exception if we are asking for information which is not present, because
     * upstroke was never detected (the voltage trace never exceeded the threshold voltage).
     */
    void CheckExceededThreshold(void);

    /**
     * Throw an exception if we are asking for information which is not present, because
     * an AP was never detected (the voltage trace never exceeded the threshold voltage
     * and then returned past it).
     */
    void CheckReturnedToThreshold(void);

public:
    /**
     * Constructor sets the data and calls CalculateProperties
     *
     * @param &rVoltage a reference to the vector of voltages
     * @param &rTime a reference to the vector of times
     * @param threshold is the threshold for determining if an AP started, defaults to -30
     *
     */
    CellProperties(const std::vector<double>& rVoltage,
                   const std::vector<double>& rTime,
                   double threshold = -30.0)
            : mrVoltage(rVoltage),
              mrTime(rTime),
              mUnfinishedActionPotentials(false),
              mThreshold(threshold)
    {
        CalculateProperties();
    }

    /**
     * Returns the maximum upstroke velocity for all APs.
     *
     * @return a vector containing the maximum upstroke velocity for all APs
     */
    std::vector<double> GetMaxUpstrokeVelocities();

    /**
     * Returns the maximum upstroke velocity for the last AP.
     * If only one incomplete AP is generated, it returns the maximal upstroke so far.
     * If the threshold is never crossed, it throws an exception.
     *
     * @return the upstroke velocity of the last AP
     */
    double GetLastMaxUpstrokeVelocity();

    /**
     * Returns the maximum upstroke velocity for the last complete AP.
     * If the threshold is never crossed, it throws an exception.
     *
     * @return the upstroke velocity of the last complete AP
     */
    double GetLastCompleteMaxUpstrokeVelocity();

    /**
     * Returns the time at which the maximum upstroke velocity occurred for all APs.
     *
     * @return a vector containing the times of maximum upstroke velocity for all APs
     */
    std::vector<double> GetTimesAtMaxUpstrokeVelocity();

    /**
     * Returns the time at which the maximum upstroke velocity for the last AP occurred.
     * If only one incomplete AP is generated, it returns the time of the maximal upstroke so far.
     * If the threshold is never crossed, it throws an exception.
     *
     * @return the time of the upstroke velocity of the last AP
     */
    double GetTimeAtLastMaxUpstrokeVelocity();

    /**
     * Returns the time at which the maximum upstroke velocity for the last complete AP occurred.
     * If the threshold is never crossed, it throws an exception.
     *
     * @return the time of the upstroke velocity of the last full AP
     */
    double GetTimeAtLastCompleteMaxUpstrokeVelocity();

    /**
     * Returns the cycle lengths for all APs.
     *
     * @return a vector containing the cycle lengths for all APs
     */
    std::vector<double> GetCycleLengths();

    /**
     * Returns the peak potentials for all APs.
     *
     * @return a vector containing the peak potentials for all APs
     */
    std::vector<double> GetPeakPotentials();

    /**
     * Returns the times of all the peak potentials for all APs.
     *
     * @return a vector containing the times of peak potentials for all APs
     */
    std::vector<double> GetTimesAtPeakPotentials();

    /**
     * Returns the time of the final peak potential.
     *
     * @return the times of the final peak potentials
     */
    double GetTimeAtLastPeakPotential();

    /**
     * Returns the last peak potential.
     *
     * @return last peak potential V_max
     */
    double GetLastPeakPotential();

    /**
     * Returns the last complete AP's peak potential.
     *
     * @return last peak potential V_max
     */
    double GetLastCompletePeakPotential();

    /**
     * Returns the time of the last complete AP's peak potential.
     *
     * @return last peak potential time.
     */
    double GetTimeAtLastCompletePeakPotential();

    /**
     * Returns the minimum resting potential BEFORE each AP.
     *
     * These are calculated as the point where the derivative
     * of the potential is lowest, i.e. when the profile
     * is flattest in between two APs.
     *
     * Or just the minimum voltage before an upstroke
     * if we can't find a flat bit.
     *
     * So this vector should be the same length as the other AP property vectors.
     *
     * @return a vector containing the resting potentials for all APs
     */
    std::vector<double> GetRestingPotentials();

    /**
     * Returns the last resting potential
     * (that would be the last entry of vector returned by #GetRestingPotentials())
     *
     * @return the resting potential BEFORE the last upstroke.
     */
    double GetLastRestingPotential();

    /**
     * Returns all the action potentials durations
     *
     * @param percentage is the repolarisation percentage that
     * the APD will be calculated for. e.g. percentage = 90 for APD90.
     *
     * @return a vector containing all the APDs
     */
    std::vector<double> GetAllActionPotentialDurations(const double percentage);

    /**
     * Returns the amplitude of the last action potential generated.
     * Throws an exception if no AP is generated.
     *
     * @param percentage is the repolarisation percentage that
     * the APD will be calculated for. e.g. percentage = 90 for APD90.
     *
     * @return the APD of the last AP
     */
    double GetLastActionPotentialDuration(const double percentage);

    /**
     * Returns the amplitude of all the action potentials calculated.
     *
     * @return a vector containing all the AP amplitudes
     */
    std::vector<double> GetActionPotentialAmplitudes();

    /**
     * Returns the amplitude of the last action potential that was recorded.
     *
     * @return The amplitude of the last action potential that was recorded.
     */
    double GetLastActionPotentialAmplitude();

    /**
     * @return a vector containing the number of above-threshold depolarisations for each Ap.
     */
    std::vector<unsigned> GetNumberOfAboveThresholdDepolarisationsForAllAps();

    /**
     * @return the number of above-threshold depolarisations for the last Ap.
     */
    unsigned GetNumberOfAboveThresholdDepolarisationsForLastAp();
};

#endif //_CELLPROPERTIES_HPP_
